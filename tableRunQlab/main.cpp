// This file is under public domain.
/*
 * Conventions:
 * Every array has 0-based index, unless calling foreign functions.
 * 
 */

/*
 * Usage: cmd [lastBG.csv energylow energyhigh]
 * -- no argument : compute from scratch
 * -- 3 argument : file name to existing background in csv, energylow, energyhigh
 *
 * The program does DFT calculation in composite fermion picture.
 * We keep the system electrically neutral all the time.
 *
 * To setup the problem, edit constants/variables in syspara.hpp.
 * To ensure successful diagonalization, please set `NsubspaceHardMax(Min)`
 * and `Energylow(high)` to appropriate values.
 *
 * Output: (all output files are prefixed with string variable /prefix/.
 * SeaFillingBG###.csv      Electron filling factor (prop to density)
 * ionBG.csv                Ion filling factor, usually also the initial density.
 * log.log                  Main log file. All of output (mainout) is copied here.
 * energy_iter.sum          Energies in each iteration recorded.
 * relDensity-diff.sum      Relative density changes in each iteration recorded.
 * ees-1                    Last set of eigenvalues calculated.
 *
 */


/*
 * The program does the following 
 * 1. Construct lattice 
 * 2. Construct green's function convolver for A_x, A_y and VCoulomb
 * 3. Load initial state
 * 4. Build hamiltonian draft, using a map
 * 5. Compile hamiltonian into some sparse matrix format.
 * 6. Iteration
 */


#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <type_traits>
#include <chrono>
#include <utility>

/*
  Notes on maps:
  - map is ordered, implemented with a balanced tree.
  - nonexistent keys has zero as initial values.
*/
#include <map>
#include <unordered_map>

#include <array>
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>//why keep it? something but not error happens during compile if deleted. but why?

using namespace std;

#include "teestream.hpp"
#include "tools.hpp"
#include "densematIO.hpp"
#include <mkl.h>
#include "commoninclude.hpp"


/* using intel MKL
   #include <mkl_types.h>
   #include <mkl.h>
*/

void premisecheck() {
	static_assert(std::is_same<MKL_INT, iint>::value, "MKL_INT != iint");  // check same type or not; print string message if false bool_constexpr
	static_assert(sizeof(MKL_Complex16)==sizeof(ComplexF64), "MKL_Complex16 != ComplexF64");
}

#include "sparsemat.hpp"
typedef MMSpmatR<iint,ComplexF64> MMSpmatR_c64;

#include "syspara.hpp"
#include "lattice.h"
#include "electromag.hpp"
#include "excorr.hpp"
#include "eigensolver.hpp"
#include "mklspmat.hpp"
//#include "DFTbuilder.hpp"

#include "extraV.hpp"
#include "Disorder.hpp"


/* main program global variables */
namespace MAIN {
	/*
	 * Size of eigenspace to allocate.
	 * First set reasonable lower bound independent of filling.
	 * Filling will be taken into account later.
	 * (Need to be) ensured to be greater than Np throughout.
	 */
	/*it's the user's responsibilty to make sure MaxSubspace is (a)large enough and (b)not exceeding machine memory to cause segmentation fault*/
	iint MaxSubspace = 6000;//min(6000, Nxy); // Will be initialized in main by setMaxSubspace. This memory will be fixed throughout once allocated; not sure why this way??? Yayun

	arma::cx_mat    Eigenvecs, EigenvecsOld;
	double			*Eigenvals=NULL;
	double          *EigenvalsOld = NULL;
	arma::mat		 SeaDensity, lastSeaDensity, newSeaDensity; // Density of all the electron cloud.
	arma::mat		 dynP_x, dynP_y, newVT, oldVT, lastVT; // V_T^* matrices and the convolvants (p-qA).
	arma::mat		 ionCoulomb;	// First used as ionDensity. Then will be set to ion potential, wiz negative of the initial VCoulomb.
	map<At,ComplexF64>	 ham_draft; // Draft space for building up hamiltonian.
	MMSpmatR_c64		 ham;

	iint Nstate_per_Elevel;
	double initElectronRadius;
	
	/*normal data record files, names of recordings*/
	const string	energyFile = prefix + "energy_iter_";
	const string	fermiFile = prefix + "fermi_";
	const string	evsdiffFilebin = prefix + "evsdiffbin_";
	const string	evsoverlapFilebin = prefix + "evsoverlapbin_";
	const string	occuAllFilebin = prefix + "occuAllbin_";
	const string	eigenAllFilebin = prefix + "eigenAllbin_";
	const string	TableOccuEigen = prefix + "TableOccuEigen_";
	const string	TableHistoryDetail = prefix + "TableHistoryDetail_";
	/*converged data of occu, ees, record files, names of recordings*/
	const string	convergedOccu = prefix + "convergedOccu";
	const string	convergedEigen = prefix + "convergedEigen";

	const string	logFile = prefix+"log.log";
	ofstream	logStream(logFile);
	//ofstream	logStream;
	//logStream.open(logFile);
	teestream	mainout(std::cout, logStream);
	//teestream   mainout(std::cout, logStream); do not know how to play with the file names
	//logStream.close();


	/*
	const string	newlogFile = prefix + "newlog.log";
	ofstream	newlogStream(newlogFile);
	teestream newmainout(std::cout, newlogStream);
	*/

	bool continueFromFile = false;
	bool continueWithVT = false;
	char *loadSeaFile, * loadVTFile = ""; // filename to read old data & continue calculation

	double* occupationNumber;
	double* occupationNumberOld;
	double fermiEnergy;

/* tests and verbose reports */
#include "main_util.hpp"
#include "DFTbuilder.hpp"

}; // namespace MAIN


/*****************************************
 * Helper functions
 *****************************************/
	

/* round up V in units of U */
inline double levelup(double V, double U) {
	return ceil(V/U)*U;
}

/* Density from Eigenvec for armadillo matrices */
void setDensity(arma::mat &rho, const arma::cx_mat &evs, const iint np) {
	using namespace arma;
	// extern iint Nx, Ny;
	arma::mat a2 = sum(pow(abs(evs.cols(0,np-1)),2),1);
	rho = reshape(a2, Nx, Ny);
}

/* Density from Eigenvec for armadillo matrices, use occupationNumber */
void setDensityTemperature(arma::mat& rho, const arma::cx_mat& evs, const int& nee_get, double* occu) {
	using namespace arma;
	// extern iint Nx, Ny;
	arma::mat a2 = pow(abs(evs.col(0)), 2)* occu[0];
	for (int i = 1;i < nee_get;i++) {
		a2 += pow(abs(evs.col(i)), 2) * occu[i];
	}	
	rho = reshape(a2, Nx, Ny);
}

/*****************************************
 * Ion distribution and initial state
 *****************************************/

/* Set a disk region of density to `filling`. `radius` is in number of grids. */
void fillDisk(arma::mat &density, const double filling, const double radius) {
	const double r2 = radius*radius;
	const double wfsq = fillingToAbs2(filling);
	const double mx = (Nx-1)/2.0, my=(Ny-1)/2.0;
	for(iint j=0; j<Ny; j++) {
		for(iint i=0; i<Nx; i++) {
			const double di=i-mx, dj=j-my;
			if(di*di+dj*dj < r2) {
				density(i,j) = wfsq;
			} else {
				density(i,j) = 0;
			}
		}
	}
}

/* Set a full Rectangular region of density to `filling`. `sideX` is in number of grids. */
void fillFullRectangle(arma::mat& density, const double filling, const int num) {
	if (num==1){//old rectangle full 
		const double sideX1 = Nx/9.0*2.5;
		const double sideX2 = Nx/9.0*6.5;
		const double QPCwidth = 10 * Nry;
		//const double GateHalfWidth = 3.0;
		const double GateHalfWidth = 1;
		//const double CornerLength = 4.0;
		const double CornerLength = 0;
		//const double GateCornerLength = 5.0;
		const double GateCornerLength = 0;
		const double sideY1 = Ny/2.0 - QPCwidth/2.0;
		const double sideY2 = Ny/2.0 + QPCwidth/2.0;
		const double wfsq = fillingToAbs2(filling);
		for (iint j = 0; j < Ny; j++) {
			for (iint i = 0; i < Nx; i++) {
				if ( abs(i-sideX1) <= GateHalfWidth*Nrx &&  j <= sideY1 || abs(i-sideX1) <= GateHalfWidth*Nrx &&  j >= sideY2 || abs(i-sideX2) <= GateHalfWidth*Nrx &&  j <= sideY1 || abs(i-sideX2) <= GateHalfWidth*Nrx &&  j >= sideY2 ) 
				{
					density(i, j) = 0;//deplete 4 gates of width GateHalfWidth*magnetic length
				}else if ( abs(i) + abs(j) <= CornerLength*(Nrx+Nry) || abs(i-Nx) + abs(j) <= CornerLength*(Nrx+Nry) || abs(i) + abs(j-Ny) <= CornerLength*(Nrx+Nry) || abs(i-Nx) + abs(j-Ny) <= CornerLength*(Nrx+Nry) ) 
				{
					density(i, j) = 0;//deplete 4 corners by CornerLength*magnetic length
				}else if ( abs(i-sideX1) + abs(j) <= CornerLength*(Nrx+Nry) || abs(i-sideX2) + abs(j) <= CornerLength*(Nrx+Nry) || abs(i-sideX1) + abs(j-Ny) <= CornerLength*(Nrx+Nry) || abs(i-sideX2) + abs(j-Ny) <= CornerLength*(Nrx+Nry) ) 
				{
					density(i, j) = 0;//deplete 8 gate corners by CornerLength*magnetic length
				}
				else {
					density(i, j) = wfsq;
				}
			}
		}
	}else if(num==2){//new period
		if (1==2){
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double yMargin = 1;
			const double QPCHalfWidth = 10;
			const double W_qpc_half = QPCHalfWidth*Nry;
			const double H = (Ny/2.0 - yMargin * Nry - W_qpc_half)/2.0;
			const double GateWidth = 1*Nry;

			for (iint i = 0; i < Nx; i++) {
				for (iint j = round( H*(cos(2*Pi/Nx*i)-1 ) + Ny/2 - W_qpc_half ); j <= round( -H*(cos(2*Pi/Nx*i)-1 ) + Ny/2 + W_qpc_half ); j++) {
						density(i, j) = wfsq;//density between the gate
				}
				for (iint j = 0; j <= round( H*(cos(2*Pi/Nx*i)-1 ) + Ny/2 - W_qpc_half - GateWidth); j++) {
						density(i, j) = wfsq;//density below the lower gate
				}
				for (iint j = Ny -1; j >= round( -H*(cos(2*Pi/Nx*i)-1 ) + Ny/2 + W_qpc_half + GateWidth); j--) {
						density(i, j) = wfsq;//density above the lower gate
				}
			}
		}else if(1==2){
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double yMargin = 1;
			const double QPCHalfWidth = 10;
			const double W_qpc_half = QPCHalfWidth*Nry;
			const double H = (Ny/2.0 - yMargin * Nry - W_qpc_half)/2.0;
			const double GateWidth = 1*Nry;

			for (iint i = 0; i < Nx; i++) {
				for (iint j = round( H*(cos(2*Pi/Nx*3*i)-1 ) + Ny/2 - W_qpc_half ); j <= round( -H*(cos(2*Pi/Nx*3*i)-1 ) + Ny/2 + W_qpc_half ); j++) {
						density(i, j) = wfsq;//density between the gate
				}
			}
		}else if(1==2){
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double yMargin = 1;
			const double QPCHalfWidth = 10;
			const double W_qpc_half = QPCHalfWidth*Nry;
			const double H = (Ny/2.0 - yMargin * Nry - W_qpc_half)/2.0;
			const double GateWidth = 1*Nry;
			const double LxPeriod = Nx/3.0;
			for (iint i = 0; i < Nx; i++) {
				for (iint j = round( H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 - W_qpc_half ); j <= round( -H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 + W_qpc_half ); j++) {
						density(i, j) = wfsq;//density between the gate
				}
				for (iint j = 0; j <= round( H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 - W_qpc_half - GateWidth); j++) {
						density(i, j) = wfsq;//density below the lower gate
				}
				for (iint j = Ny -1; j >= round( -H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 + W_qpc_half + GateWidth); j--) {
						density(i, j) = wfsq;//density above the lower gate
				}
			}	
		}else if(1==2){
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double yMargin = 0.5;
			const double QPCHalfWidth = 5;
			const double W_qpc_half = QPCHalfWidth*Nry;
			const double H = (Ny/2.0 - yMargin * Nry - W_qpc_half)/2.0;
			const double GateWidth = 1*Nry;
			const double LxPeriod = Nx/3.0;
			for (iint i = round(LxPeriod/2.0); i < Nx - round(LxPeriod/2.0); i++) {
				for (iint j = round( H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 - W_qpc_half ); j <= round( -H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 + W_qpc_half ); j++) {
						density(i, j) = wfsq;//density between the gate
				}
				for (iint j = 0; j <= round( H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 - W_qpc_half - GateWidth); j++) {
						density(i, j) = wfsq;//density below the lower gate
				}
				for (iint j = Ny -1; j >= round( -H*(cos(2*Pi/LxPeriod*i)-1 ) + Ny/2 + W_qpc_half + GateWidth); j--) {
						density(i, j) = wfsq;//density above the lower gate
				}
			}
			for (iint i = 0; i < round(LxPeriod/2.0); i++) {
				for (iint j = 0; j <= Ny -1 ; j++) {
						density(i, j) = wfsq;//density in left lead
				}
			}
			for (iint i = Nx - round(LxPeriod/2.0); i <= Nx - 1; i++) {
				for (iint j = 0; j <= Ny -1 ; j++) {
						density(i, j) = wfsq;//density in right lead
				}				
			}	
		}else if(1==2){
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double torelance = 0.5;//in units of 1
			const double QPCWidth = 10.0;
			const double W_qpc = QPCWidth * Nry;
			const double H = (Ny - W_qpc)/2.0;
			const double GateWidth = 3.0 * Nry;
			const double leftPos = round(Nx * 7.0/24.0);//yayun,20240108, need a round here. otherwise H can not be reached when gatewidth too small and leftpos is non-integral
			const double rightPos = Nx - leftPos;
			for (iint i = 0; i < Nx ; i++) {
				for (iint j = round( H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j <= round( Ny-1 - H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j++) {
						density(i, j) = wfsq;//density between the gate
				}				
			}	
		}else if(1==2){
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double torelance = 0.1;//in units of 1
			const double QPCWidth = 6;
			const double W_qpc = QPCWidth * Nry;
			const double H = (Ny - W_qpc)/2.0;
			const double GateWidth = 3.0 * Nry;
			const double leftPos = round(Nx * 7.0/24.0);//yayun,20240108, need a round here. otherwise H can not be reached when gatewidth too small and leftpos is non-integral
			const double rightPos = Nx - leftPos;
			//cout << "H: " << H << endl;
			//cout << "W_qpc: " << W_qpc << endl;
			for (iint i = 0; i < Nx ; i++) {
				for (iint j = round( H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j <= round( Ny-1 - H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j++) {
						density(i, j) = wfsq;//density between the gate
				}				
			}	
		}else if(1==1){//To do something to make sure particle = degeneracy ,yayun,240118
			const double wfsq = fillingToAbs2(filling);
			density.zeros();
			const double torelance = 0.1;//in units of 1
			const double QPCWidth = 20;
			const double W_qpc = QPCWidth * Nry;
			const double H = (Ny - W_qpc)/2.0;
			const double GateWidth = 3.0 * Nry;
			const double leftPos = round(Nx * 7.0/24.0);//yayun,20240108, need a round here. otherwise H can not be reached when gatewidth too small and leftpos is non-integral
			const double rightPos = Nx - leftPos;
			const int edgeMargin = 1.0 * Nry;
			//cout << "H: " << H << endl;
			//cout << "W_qpc: " << W_qpc << endl;
			for (iint i = 0; i < Nx ; i++) {
				for (iint j = round( H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j <= round( Ny-1 - H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j++) {
						density(i, j) = wfsq;//density between the gate
				}				
			}
			const double marginDensity = 0.00000001;
			//left and right edge margin
			for (iint j = 0; j < Ny ; j++) {
				for (iint i = 0; i <= edgeMargin-1 ; i++) {				
					density(i, j) = marginDensity;//density left edge margin
				}
				for (iint i = Nx - edgeMargin; i <= Nx-1 ; i++) {
					density(i, j) = marginDensity;//density right edge margin
				}
			}

			//up and down gate margin
			for (iint i = 0; i < Nx ; i++) {
				for (iint j = round( H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j <= round( H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )) + edgeMargin-1; j++) {
						density(i, j) = marginDensity;//down gate
				}
				for (iint j = round( Ny-1 - H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) ))+1 - edgeMargin; j <= round( Ny-1 - H * (exp(-0.5*pow((i-leftPos)/GateWidth,2)) + exp(-0.5*(pow((i-rightPos)/GateWidth,2))) )); j++) {
						density(i, j) = marginDensity;
				}				
			}	
		}
	}
}


/* add depleting region barrier */
void addDeltaBarrierAndGateBarrier(arma::mat& DeltaBarrierPotential, const enum SHAPE shapeCS, const double DeltaBarrierNp, const double geometryYoverXDeltaBarrier) {
	DeltaBarrierPotential.zeros();
	const double DeltaBarrierHeight = 1000.0;

	if (shapeCS == Circle) {
	   
	}
	else if (shapeCS == Rectangle) {
		auto xLenDeltaBarrier = sqrt(DeltaBarrierNp / fillingToAbs2(ionFilling) / geometryYoverXDeltaBarrier); // length of side x region for the initial state and the ion background.
		auto yLenDeltaBarrier = xLenDeltaBarrier * geometryYoverXDeltaBarrier;
		auto yLenDeltaBarrierOpen = yLenDeltaBarrier * 1 / 4;
		double DeltaBarrierWidth = 0.0;// l_{B_{ reference }}
		int roundWidth = 0;//width as grid number
		//find the values of grid edge
	   const int xgridvaluemin = round((Nx - xLenDeltaBarrier) / 2.0);
	   const int xgridvaluemax = round((Nx + xLenDeltaBarrier) / 2.0);
	   const int ygridvaluemin = round((Ny - yLenDeltaBarrier) / 2.0);
	   const int ygridvaluemax = round((Ny + yLenDeltaBarrier) / 2.0);
	   const int ygridvalueOpenmin = round((Ny - yLenDeltaBarrierOpen) / 2.0);
	   const int ygridvalueOpenmax = round((Ny + yLenDeltaBarrierOpen) / 2.0);
	   for (iint i = xgridvaluemin; i <= xgridvaluemax; i++) {
		   roundWidth = round(Nrx * DeltaBarrierWidth) - 1;
		   if (roundWidth < 0) { roundWidth = 0; }
		   for(iint iw=0; iw<= roundWidth ; iw++){
		       DeltaBarrierPotential(i, ygridvaluemin - iw) = DeltaBarrierHeight;
		       DeltaBarrierPotential(i, ygridvaluemax + iw) = DeltaBarrierHeight;
		   }
	   }
	   for (iint j = ygridvaluemin + 1; j <= ygridvalueOpenmin; j++) {
		   roundWidth = round(Nry * DeltaBarrierWidth) - 1;
		   if (roundWidth < 0) { roundWidth = 0; }
		   for (iint iw = 0; iw <= roundWidth; iw++) {
			   DeltaBarrierPotential(xgridvaluemin - iw, j) = DeltaBarrierHeight;
			   DeltaBarrierPotential(xgridvaluemax + iw, j) = DeltaBarrierHeight;
		   }
	   }
	   for (iint j = ygridvalueOpenmax; j <= ygridvaluemax-1; j++) {
		   roundWidth = round(Nry * DeltaBarrierWidth) - 1;
		   if (roundWidth < 0) { roundWidth = 0; }
		   for (iint iw = 0; iw <= roundWidth; iw++) {
			   DeltaBarrierPotential(xgridvaluemin - iw, j) = DeltaBarrierHeight;
			   DeltaBarrierPotential(xgridvaluemax + iw, j) = DeltaBarrierHeight;
		   }
	   }
	}
	else if(shapeCS == FullRectangle ){
		fillFullRectangle(DeltaBarrierPotential, 1/*any nonzero number*/, 2/*means using old rectangle full or new period*/);
		for (iint j = 0; j < Ny; j++) {
			for (iint i = 0; i < Nx; i++) {
				if (DeltaBarrierPotential(i,j) == 0){
					DeltaBarrierPotential(i,j) = DeltaBarrierHeight;
				}else{
					DeltaBarrierPotential(i,j) = 0;
				}
			}
		}
		const int mx = (Nx - 1) / 2.0, my = (Ny - 1) / 2.0;//set center position
		const int W_imp = 2 * Nry;
		for (iint j = my - W_imp; j <= my + W_imp; j++) {
			for (iint i = mx - W_imp; i <= mx + W_imp; i++) {
			//	if (DeltaBarrierPotential(i,j) == 0){
					DeltaBarrierPotential(i,j) = 0.0;
				//}else{
				//	DeltaBarrierPotential(i,j) = 0;
			//	}
			}
		}

	}

}


/* Set a Rectangular region of density to `filling`. `sideX` is in number of grids. */
void fillRectangle(arma::mat& density, const double filling, double sideX, double sideY) {
	sideX = sideX/2.0; 
	sideY = sideY/2.0;
	const double wfsq = fillingToAbs2(filling);
	const double mx = (Nx - 1) / 2.0, my = (Ny - 1) / 2.0;//set center position
	for (iint j = 0; j < Ny; j++) {
		for (iint i = 0; i < Nx; i++) {
			const double di = i - mx, dj = j - my;
			if ( abs(di)<= sideX &&  abs(dj)<= sideY) {
				density(i, j) = wfsq;
			}
			else {
				density(i, j) = 0;
			}
		}
	}
}



/*
 * Prologue prints the parameters before the serious computaions are done.
 */
void prologue() {
	using namespace MAIN;
	using namespace std;

	mainout<<"========================================"<<endl
	    <<"DFT calculation begin"<<endl
	       <<"========================================"<<endl<<endl;

	
	if(continueFromFile) {
		mainout << "Continuing from density profile: "<<loadSeaFile<<endl<<endl;
		mainout << "Continuing from VT profile: " << loadVTFile << endl
			<< endl;
	}
	

	mainout<<"System size in plain magnetic length: "<<Lx<<" by "<<Ly<<endl
	       <<"  on a rectangular sample." << endl
	       <<"Magnetic length resolution is "<<Nrx<<" by "<<Nry<<endl
	       <<"Total discretization grid is "<<Nx<<"x"<<Ny<<"="<<Nx*Ny<<endl
	       <<"Filling to wave function density factor (for delta B =0) is "<<NVSTRING(filling_vs_abs2)<<endl
	       <<endl;

	mainout << "Given ion filling: " << ionFilling << ',' << endl;
	       //<<"The resolution per reduced magnetic length is "<<Nrxu<<"x"<<Nryu<<'.'<<endl<<endl;

	mainout<<"Given the following parameters:"<<endl
	       <<'\t'<<NVSTRING(alpha_m)<<endl
	       <<'\t'<<NVSTRING(Dielec)<<endl
	       <<'\t'<<NVSTRING(KE_scale)<<endl
	       <<"The cyclotron gap is predicted to be: "<<CyclotronGap<<'.'<<endl<<endl;

	mainout<<"Anisotropy parameters: "<<endl
	       <<'\t'<<NVSTRING(Mx_vs_My)<<endl
	       <<'\t'<<NVSTRING(DielecX_vs_DielecY)<<endl<<endl;

	mainout<<"Hamiltonian switches: "<<endl
	       <<'\t'<<NVSTRING(scaleHartree)<<endl
	       <<'\t'<<NVSTRING(scaleIon)<<endl
	       <<'\t'<<NVSTRING(scaleXC)<<endl
	       <<'\t'<<NVSTRING(scaleVT)<<endl;

	mainout<<"Ions fill the central disk of grid radius "<<ionRadius<<'.'<<endl
	       <<"The system has "<<Np<<" particles."<<endl
	       <<"Maximum eigensubspace size "<<MaxSubspace<<endl<<endl;

	if(!continueFromFile){
		mainout<<"Electons has an initial filling "<<initElectronFilling<<','<<endl
		       <<"with radius "<<initElectronRadius<<'.'<<endl
		       <<endl;
	}
	
	mainout <<"Confinement "<<confinementName<<" is added; height is "<<confinementHeight<<'.'<<endl;

	mainout<<"Predicted number of states per lambda level is "<<Nstate_per_Elevel<<'.'<<endl;
	mainout<<"(Without confinement, across the whole grid, it is "<<Nstate_per_lambda<<".)"<<endl;

	mainout<<"Energy bound guesses: " << Energylow <<" -- "<<Energyhigh << endl;

	mainout<<endl<<"========================================"<<endl<<endl;
}


int main(int argc, char** argv) {

	using namespace TicToc;
	using namespace COULOMB;
	using namespace VECTORPOT;
	using namespace VTSTAR;
	using namespace MAIN;


	premisecheck();
	auto startTime = high_resolution_clock::now();
	//ofstream densityDiffout(densityDiffFile);
	//ofstream densityDiffoutAim(densityDiffFileAim);
	//ofstream VTDiffout(VTDiffFile);
	//ofstream VTDiffoutAim(VTDiffFileAim);


	/*****************************************
	 * Parse command line arguments. See file header for usage.
	 *****************************************/
	argc--;// function name "main" is also counted as a command line argument
	switch (argc) {//argc is the number of input parameters
	case 0:
		continueFromFile = false;
		break;
	case 3:
		continueFromFile = true;
		loadSeaFile = argv[1];
		Energylow = atof(argv[2]);
		Energyhigh = atof(argv[3]);
		break;
	case 4:
		//parse_config(argv[1]);
		continueFromFile = true;
		loadSeaFile = argv[1];
		loadVTFile = argv[2];
		continueWithVT = true;
		Energylow = atof(argv[3]);
		Energyhigh = atof(argv[4]);
		//manualEnergyBound = true;
		break;
	default:
		std::cerr << "Wrong usage. " + string(argv[0]) + " [lastBG.csv energylow energyhigh]" << endl;
		exit(-1);
	}

	/*****************************************
	 * System allocations and initializations.
	 * -- Allocate matrices within namespace MAIN.
	 * -- Read old filling density file if applicable. Set up initial ion & electron density.
	 * -- Estimate parameters for efficient computation of eigensystem.
	 * ----  Estimate Nstate_per_Elevel, heuristically estimate MaxSubspace, Energyhigh
	 * -- Print prologue i.e. system info.
	 * -- Compute Greens functions i.e. convolution kernels.
	 * -- Construct initial hamiltonian.
	 * -- Compute initial potential energy value, i.e. without V_T*, kinetic term
	 *****************************************/

	 /*make iteration parameters*/
	int countRho = 0;
	SeaDensity = arma::mat(Nx, Ny, arma::fill::zeros);
	lastSeaDensity = arma::mat(Nx, Ny, arma::fill::zeros);
	newSeaDensity = arma::mat(Nx, Ny, arma::fill::zeros);
	ionCoulomb = arma::mat(Nx, Ny, arma::fill::zeros);
	dynP_x = arma::mat(Nx, Ny, arma::fill::zeros);
	dynP_y = arma::mat(Nx, Ny, arma::fill::zeros);
	newVT = arma::mat(Nx, Ny, arma::fill::zeros);
	oldVT = arma::mat(Nx, Ny, arma::fill::zeros);
	lastVT = arma::mat(Nx, Ny, arma::fill::zeros);
	double Etotal = 0.0, KE = 0.0, Eext = 0.0, EHartree = 0.0, Excsum = 0.0;
	double relDensityDiff = 0.0, relDensityDiffAim = 0.0, relVTDiff = 0.0, relVTDiffAim = 0.0;
	/*sum up total potential to estimate energy bounds*/
	arma::mat totalPotentialVxc(Nx, Ny, arma::fill::zeros);
	/* initialize FailureAndBreak, to do/add */
	arma::umat TableOrder, TableOrdertempInverse, TableOrderUseOldResults;
	MakeTableOrderDimension(TemperatureNum, MySecondSetNum, TableOrder, TableOrdertempInverse, TableOrderUseOldResults);
	//mainout << TemperatureNum << "   " << MySecondSetNum <<"  " <<endl ;
	TableOrder.save(prefix + "TableOrder.csv", arma::csv_ascii);
	/*read or make files for table parameters*/
	arma::mat TableAllinOne(TemperatureNum*MySecondSetNum+1, 35, arma::fill::ones);
	TableAllinOne = TableAllinOne * (-1.0);//initialize to -1;
	//update tables: its size
	TableAllinOne(TemperatureNum * MySecondSetNum, 4) = TemperatureNum;
	TableAllinOne(TemperatureNum * MySecondSetNum, 5) = MySecondSetNum;
	TableAllinOne(TemperatureNum * MySecondSetNum, 6) = save_cycleNum1;//density
	TableAllinOne(TemperatureNum * MySecondSetNum, 7) = save_cycleNum2;//potential
	//TableAllinOne.save(prefix + "TableAllinOneA.csv", arma::csv_ascii);
	//arma::mat TableOccuEigen(2, OccuEigenStorage, arma::fill::zeros);//static_cast<int>(TemperatureNum * MySecondSetNum * (int)2)
	//arma::mat TableHistoryDetail(4, OccuEigenStorage, arma::fill::zeros);	

	int LetsStartWhereStops = 0;// the key idea is to only decide which LetsStartWhereStops to use and check file existence(to be done); which round to use is decided in for loop 
	const string alltable = prefix + "TableAllinOne.csv";
	if (KeepUsingInitialConditions == 1) {
		if (exists_test(alltable) == false) {//TableAllinOne not found
			mainout << "Couldn't read file: TableAllinOne.csv" << endl;
			if (exists_test(TableOccuEigen) == true || exists_test(TableHistoryDetail) == true) {//but others found
				mainout << "You might forget the TableAllinOne.csv file, aborting" << endl;
				exit(-2);
			}
			else
			{
				mainout << "No existing data to proceed from, start from scratch." << endl;
				LetsStartWhereStops = 0;//others not found either; Start from zero of table para 
			}
		}
		else if (!TableAllinOne.load(alltable, arma::auto_detect)) {//TableAllinOne found but can not load
			mainout << "Fail to load " << alltable << ", aborting." << endl;
			exit(-2);
		}else{//TableAllinOne found and loaded
			mainout << "TableAllinOne found and loaded. " << endl;
			if (1==2 /*do a full file check*/ /*TableAllinOne(TemperatureNum * MySecondSetNum, 0) >= 1 && (exists_test(TableOccuEigen) == false || exists_test(TableHistoryDetail) == false)*/) 
			{
				mainout << "You might forget some bin files, aborting" << endl;
				exit(-2);
			}
			else{
				//TableAllinOne.save(prefix + "TableAllinOneB.csv", arma::csv_ascii);
				int tableNumFinished = TableAllinOne(TemperatureNum * MySecondSetNum, 0);//read LetsStartWhereStops
				int tableNumStarted = TableAllinOne(TemperatureNum * MySecondSetNum, 1);//read LetsStartWhereStops
				if (tableNumStarted == -1) {//not even the first step started
					LetsStartWhereStops = 0;//Start from zero of table para
				}
				else {
					if (tableNumStarted == tableNumFinished) {//previous calculation just finished tableNumStarted and stops right there before starting the next
							LetsStartWhereStops = tableNumFinished + 1;//start the next calculation
					}else if (tableNumStarted == tableNumFinished + 1) {//previous calculation in the middle of iterations; caution: tricky and risky condition, -1+1=0				   
						LetsStartWhereStops = tableNumStarted;
						/*search for storage for initial conditions of iterations will be done in the for loop*/
					}
					else {
						mainout << "not a reasonable tableNumStarted or tableNumFinished, aborting" << endl;
						exit(-2);
					}
				}
				mainout << "tableNumFinished: " << tableNumFinished << ", tableNumStarted:" << tableNumStarted << endl;
			}
			if (LoadAndReiterateConvergedDada == 0) {
				mainout << "should KeepUsingInitialConditions=1 with LetsStartWhereStops, " << LetsStartWhereStops << " , and use it because LoadAndReiterateConvergedDada=0" << endl;
			}
			else if(LoadAndReiterateConvergedDada == 1) {
				mainout << "should KeepUsingInitialConditions=1 with LetsStartWhereStops, " << LetsStartWhereStops << " , but use LetsStartWhereStops = 0 because LoadAndReiterateConvergedDada=1" << endl;
				LetsStartWhereStops = 0;
			}
		}		
}
	mainout << "LetsStartWhereStops: " << LetsStartWhereStops << endl;
	/*
	 * Read old filling density file if applicable.
	 * Set up initial ion & electron density.
	 */
	if (continueFromFile) {
		if (!SeaDensity.load(loadSeaFile, arma::auto_detect)) {
			mainout << "Couldn't read density profile: " << loadSeaFile << endl
				<< "Aborting." << endl;
			exit(-2);
		}
		if (continueWithVT)
		{
			if (!oldVT.load(loadVTFile, arma::auto_detect))
			{
				mainout << "Couldn't read VT profile: " << loadVTFile << endl
					<< "Continue without." << endl;
				oldVT.fill(0.0); // caution against partial overwrite
				continueWithVT = false;
			}
			else {
				oldVT.save(prefix + "VT0con.csv", arma::csv_ascii);//write VT0con.csv only when it is continued
				oldVT.save(prefix + "VT0.csv", arma::csv_ascii);//write both cause later will skip writing this VT0
				lastVT = oldVT;
			}
		}
		SeaDensity.save(prefix + "SeaFillingBG0con.csv", arma::csv_ascii);//write SeaFillingBG0con.csv only when it is continued
		SeaDensity = fillingToAbs2(SeaDensity);


		/* ensure ion density has integer number of electrons */
		Np = round(accu(SeaDensity));

		if (shape == Circle) {
			ionRadius = sqrt(Np / fillingToAbs2(ionFilling) / Pi);
			if (2 * ionRadius > min(Nx, Ny)) {
				mainout << "Lattice array size not big enough for the ion disk." << endl
					<< NVSTRING(ionRadius) << endl
					<< NVSTRING(Nx) << endl
					<< NVSTRING(Ny) << endl
					<< NVSTRING(Np) << endl
					<< "Aborting." << endl;
				exit(1);
			}
			fillDisk(ionCoulomb, ionFilling, ionRadius);
		}
		else if (shape == Rectangle) {
			xLen = sqrt(Np / fillingToAbs2(ionFilling) / geometryYoverX);
			yLen = xLen * geometryYoverX;
			if (xLen > Nx || yLen > Ny) {
				mainout << "Lattice array size not big enough for the ion Rectangle." << endl
					<< NVSTRING(xLen) << endl
					<< NVSTRING(yLen) << endl
					<< NVSTRING(Nx) << endl
					<< NVSTRING(Ny) << endl
					<< NVSTRING(Np) << endl
					<< "Aborting." << endl;
				exit(1);
			}
			fillRectangle(ionCoulomb, ionFilling, xLen, yLen);
		}
		else if (shape == FullRectangle) {
			fillFullRectangle(ionCoulomb, ionFilling, shape_fill_num);
			Np = accu(ionCoulomb);//accumulate to have real electron number, is int
		}
		//double tNp = accu(ionCoulomb);//accumulate
		//ionCoulomb.transform([=](double x) {return x * Np / tNp; });//ionCoulomb not integer is ok in the following

	}
	else { // not continue with existing files, but set up densities of ion and electron from parameters
	 // here ionCoulomb is the ion density distribution
		if (shape == Circle) {
			if (2 * ionRadius > min(Nx, Ny)) {
				mainout << "Lattice array size not big enough for the ion disk." << endl
					<< NVSTRING(ionRadius) << endl
					<< NVSTRING(Nx) << endl
					<< NVSTRING(Ny) << endl
					<< NVSTRING(Np) << endl
					<< "Aborting." << endl;
				exit(1);
			}
			fillDisk(ionCoulomb, ionFilling, ionRadius);
		}
		else if (shape == Rectangle) {
			if (xLen > Nx || yLen > Ny) {
				mainout << "Lattice array size not big enough for the ion Rectangle." << endl
					<< NVSTRING(xLen) << endl
					<< NVSTRING(yLen) << endl
					<< NVSTRING(Nx) << endl
					<< NVSTRING(Ny) << endl
					<< NVSTRING(Np) << endl
					<< "Aborting." << endl;
				exit(1);
			}
			fillRectangle(ionCoulomb, ionFilling, xLen, yLen);
		}
		else if (shape == FullRectangle) {
			fillFullRectangle(ionCoulomb, ionFilling, shape_fill_num);
			//ionCoulomb.fill(fillingToAbs2(ionFilling));//uniform back ground filling, which might help attract density close to gates?
			Np = accu(ionCoulomb);//accumulate to have real electron number
			mainout << "-- tp and Np." << Np << endl;
		}
		fillFullRectangle(SeaDensity, ionFilling * 1.0, shape_fill_num);//initialize the electron density to perfect screening
		
		//SeaDensity.fill( fillingToAbs2(ionFilling * 0.9 ) );

		/* the number of electrons is integer, which is ensured by definition of Np as int type*/
		double tNp = accu(ionCoulomb);
		//SeaDensity.transform([=](double x) {return x * Np / tNp; });//use correct number density as initial density
		//ionCoulomb.transform([=](double x) {return x * Np / tNp; });
		mainout << "-- tNp : " << tNp << endl;
		Np = accu(SeaDensity);
		mainout << "-- Np : " << Np << endl;
		/*
		if (1==2){
			if (shape == Circle) {
				initElectronRadius = ionRadius * sqrt(ionFilling / initElectronFilling);
				if (2 * initElectronRadius > min(Nx, Ny)) {
					mainout << "Lattice array size not big enough for the electron disk." << endl
						<< NVSTRING(initElectronRadius) << endl
						<< NVSTRING(Nx) << endl
						<< NVSTRING(Ny) << endl
						<< NVSTRING(Np) << endl
						<< "Aborting." << endl;
					exit(1);
				}
				fillDisk(SeaDensity, initElectronFilling, initElectronRadius);//  0110 test density
				//fillDisk(SeaDensity, 0.0, initElectronRadius);
			}
			else if (shape == Rectangle) {
				if (xLen > Nx || yLen > Ny) {
					mainout << "Lattice array size not big enough for the electron Rectangle." << endl
						<< NVSTRING(xLen) << endl
						<< NVSTRING(yLen) << endl
						<< NVSTRING(Nx) << endl
						<< NVSTRING(Ny) << endl
						<< NVSTRING(Np) << endl
						<< "Aborting." << endl;
					exit(1);
				}
				fillRectangle(SeaDensity, initElectronFilling, xLen, yLen);
			}
			else if (shape == FullRectangle) {
				fillFullRectangle(SeaDensity, initElectronFilling, shape_fill_num);
				double Np = accu(SeaDensity);//accumulate to have real electron number. 231226, Yayun, why we have this test anymore?
			}
		}
		*/
	}		// here SeaDensity is the electron density distribution; no further processing because it will be determined by orbitals later
	abs2ToFilling(ionCoulomb).save(prefix + "ionBG.csv", arma::csv_ascii);
	abs2ToFilling(SeaDensity).save(prefix + "SeaFillingBG0.csv", arma::csv_ascii);//This might be different from SeaFillingBG0con.csv

	/*
	 * Estimate parameters for efficient computation of eigensystem.
	 */

	 /* Estimate Nstate per lambda level over the ion region (disk/rectangle) only */
	/*
	if (shape == Circle) {
		Nstate_per_Elevel = Nstate_per_lambda * Pi * pow(ionRadius, 2) / 1.0 / Nx / Ny;
	}
	else if (shape == Rectangle) {
		Nstate_per_Elevel = Nstate_per_lambda * Pi * xLen * yLen / 1.0 / Nx / Ny;
	}

	if (Nstate_per_Elevel == 0) { // when, say, box level spacing is larger than landau level spacing. Use a dummy value.
		Nstate_per_Elevel = 1;
	}
	else if (Nstate_per_Elevel > Nstate_per_lambda) { // when ionRadius > Nx, etc
		Nstate_per_Elevel = Nstate_per_lambda;
	}
	*/

	prologue();
	/*
	 * Compute Greens functions i.e. convolution kernels.
	 */
	initGreenA();
	createConvTaskA();
	initGreenCoulomb();//DielecX_vs_DielecY is considered when creating green function 1/r
	createConvTaskCoulomb();
	initGreenVT();
	createConvTaskVT();
	mainout << "-- A & Vcoulomb & VT green's function computed." << endl;
	int Nsubspace = min(MaxSubspace, 3000);//feast parameter: estimated number of eigenstates in the range
	int EfMode = 1; /*0:use all orbitals (which is Nee) ; 1:use up to the first orbital that is occupied less than C_min (which is fixNum); 2: fix Ef;
					here operations can be crude, real fix later; */
	int fixNum = 0; //if EfMode=1, fixNum is the orbital #;if EfMode=2, fixNum=nee_get is the orbital # 
	double lastEmin = 0.0;//the lowest eigen energy
	double nee_wantdecimal = 0.0;//sum of occupation numbers
	/* start anealling   */
	for (auto tableNum = LetsStartWhereStops; tableNum <= static_cast<int>(TableOrder.n_rows) - 1; tableNum++) {
		if (TableAllinOne(tableNum, 19) == 1) {
			if (LoadAndReiterateConvergedDada == 0) {
				mainout << "tableNum " << tableNum << " status," << TableAllinOne(tableNum, 19) << " ,next" << endl;
				continue;//this means this tableNum already converged; Converged = 1; has been started = 0; finished but not converged = -2; default = -1;
			}
			else if (LoadAndReiterateConvergedDada == 1) {
				mainout << "tableNum " << tableNum << " status," << TableAllinOne(tableNum, 19) << " ,reiterate" << endl;
			}
		}
		bool continueWithVTfromTable = false;//this will be later changed if use existing files
		//load table parameters
		int ikBT = TableOrder(tableNum, 0);
		int id_set = TableOrder(tableNum, 1);
		kBT = MyTemperatureSet(ikBT);
		BExternalGlobal = MySecondSet(id_set);
		//KE_scale = KE_scaleReferenceB / sqrt(1.0 + BExternalGlobal);//0502Yayun consider change of CF mass with the external magnetic field;
		//Yayun, 240117, let's make sure other factors not change, and focus on phase
		mainout << "ikBT: " << ikBT << ", kBT: " << kBT << ", id_set: " << id_set << ", value: " << MySecondSet(id_set) << endl;
		//load initial conditions
		
		//initial update tables
		mainout << "TableAllinOne.size(): " << TableAllinOne.size() << endl;
		/*there might be a overwritten by the same values , each time tableNum is entered here*/
		TableAllinOne(TemperatureNum* MySecondSetNum, 1) = tableNum;/*tableNumStarted, this is the tableNum that has been started; this will be saved only if iter>=1 entered,
        otherwise no effect; so will not return to previous stage once saved;*/
		TableAllinOne(tableNum, 23)= MyTemperatureSet(ikBT);
		TableAllinOne(tableNum, 24) = MySecondSet(id_set);
		TableAllinOne(tableNum, 25) = Lx;TableAllinOne(tableNum, 26) = Ly;TableAllinOne(tableNum, 27) = Nrx;TableAllinOne(tableNum, 28) = Nry;
		//read from tables
		int roundNum=TableAllinOne(tableNum, 14);//this is recorded if this tableNum has been started previously
		int iterLast = TableAllinOne(tableNum, 10);//this is recorded if this tableNum has been started previously
		std::string TableRoundNum = "t" + std::to_string(tableNum) + "r" + std::to_string(roundNum);//old one
		mainout << "tableNum: "<< tableNum << ", roundNum: " << roundNum << ", iterLast: " << iterLast << endl;
		if (tableNum >= 1 /*must have saved data of previous tableNum<=0*/ || roundNum >= 1 /*must have saved data of this or previous roundNum*/
			|| (roundNum ==0 && iterLast>=1) /*make sure this roundNum ==0 contains useful data*/ ) {
			if (roundNum == -1 || (roundNum == 0 && (iterLast == -1 || iterLast ==0) )) {//roundNum==-1 means never touched this tableNum before
				if (loadPreviousTableResults(tableNum, SeaDensity, oldVT, Energylow, Energyhigh, Nsubspace, TableOrdertempInverse, TableOrderUseOldResults, 
					TableAllinOne, lastEmin)) {
					roundNum += 1;
					mainout<< "use previous table" <<endl;
					continueWithVTfromTable = true;
				}
				else {
					mainout << "data of previous tableNum not found. Aborting. " << endl;
					exit(1);
				}
			}
			else{//there should be some data in previous/present round
				if (loadPreviousRoundResults(tableNum, SeaDensity, oldVT, Energylow, Energyhigh, Nsubspace, TableOrder, TableAllinOne, lastEmin)==0) {
					//mainout << "Energylow: " << Energylow << ", Energyhigh: " << Energyhigh << endl;
					roundNum += 1;//this way roundNum is taken care of automatically, i.e., increase by 1 
					mainout << "use previous round" << endl;
					continueWithVTfromTable = true;
				}
				else {
					mainout << "data of previous round not found. Aborting. " << endl;
					exit(1);
				}				
			}
		}
		else{ 
			roundNum += 1; }//otherwise will continue with nothing but input arguments or pre-set initial conditions
		//TableAllinOne(ikBT, id_set, 14) = roundNum;//not update here because nothing significant happen yet
		//TableAllinOne(tableNum, 10) = -2;//set to -2 to indicate that it is at this step, no real iteration yet but started; not saved, and changed when saved, so meaningless
		TableAllinOne(tableNum, 14) = roundNum;// the zeroth step finishes; it will still be unchanged as before if iter>=1 is not entered, otherwise will be recorded
		TableAllinOne(tableNum, 19) = 0;/* the zeroth step finishes; it will still be unchanged as before if iter>=1 is not entered, otherwise will be recorded
		Converged = 1; has been started = 0; finished but not converged = -2; default = -1; */
		TableRoundNum = "t" + std::to_string(tableNum) + "r" + std::to_string(roundNum);//new one
		mainout << "TableRoundNum: "<< TableRoundNum << endl;
		//these are arrays of recent historical data considered for convergence
		arma::vec AllbigDiff(lookBackStepNum, arma::fill::zeros);
		arma::vec AlldensityDiff(lookBackStepNum, arma::fill::zeros);
		arma::vec AllVTDiff(lookBackStepNum, arma::fill::zeros);
		arma::vec AllkineticEnergy(lookBackStepNum, arma::fill::zeros);
		arma::vec Allnewparticle(lookBackStepNum, arma::fill::zeros);
		/*relocate arrays to make sure it contains the new parameters before tailor their size again later on*/
		if (tableNum == LetsStartWhereStops){
			/*yayun 20240219, segmentation fault repeated occurs for 2nd zeroth diagonalization,
			see E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\80
			here try use useGuessSubspaceZeroth = true to test if helpful. So not erase Eigenvec. */
			Eigenvecs = arma::cx_mat(Nxy, MaxSubspace, arma::fill::zeros);
			Eigenvals = new double[MaxSubspace];// 
			occupationNumber = new double[MaxSubspace];
		}
		mainout << "roundNum: " << roundNum << endl;
		if (externalPotentialGivenFromOutside == false) {//this means to make external potentials in the following code
			if (tableNum == 0 && roundNum==0 && TableAllinOne(tableNum, 10) == -1 /*This indicates purely fresh calculations, onhy then need update potentials*/) {
				//add disorder
				arma::mat disorder(Nx, Ny, arma::fill::zeros);
				arma::mat DeltaBarrierPotential(Nx, Ny, arma::fill::zeros);
				uniformdis(disorder);
				//longRangeImpurity(disorder);
				calcCoulomb(ionCoulomb);
				if (useAnalyticBackgroundPotential == true) {
					if (shape == Circle) {
						fillDisk(ionCoulomb, ionFilling, ionRadius);
					}
					else if (shape == Rectangle) {
						fillRectangle(ionCoulomb, ionFilling, xLen, yLen);
					}
				}
				ionCoulomb = -COULOMB::VCoulomb * scaleIon; //Now ionCoulomb is the ion potential instead, DielecX_vs_DielecY already included in green function 1/r definition
				/*External potential treated seperately*/
				/*add delta potential, merge into ionCoulomb*/
				addDeltaBarrierAndGateBarrier(DeltaBarrierPotential, shape, insideDeltaBarrierNp, geometryYoverX);
				ionCoulomb.save(prefix + "ionCoulombAlways.csv", arma::csv_ascii);//saved is used
				disorder.save(prefix + "disorder.csv", arma::csv_ascii);
				DeltaBarrierPotential.save(prefix + "DeltaBarrier.csv", arma::csv_ascii);
				mainout << "total Coulomb ionCoulomb: " << accu(ionCoulomb)<< ", total disorder"<< accu(disorder) << endl;
				ionCoulomb = ionCoulomb + disorder * disorderStrength + DeltaBarrierPotential;
				disorder.clear();//delete disorder because not used anymore 
				DeltaBarrierPotential.clear();
				mainout << "disorder and external potential added. " << endl;
			}
			else{//not purely fresh calculations
				arma::mat disorder(Nx, Ny, arma::fill::zeros); arma::mat DeltaBarrierPotential(Nx, Ny, arma::fill::zeros);
				const string disorderString = prefix + "disorder.csv", ionCoulombString = prefix + "ionCoulombAlways.csv", 
					DeltaBarrierString = prefix + "DeltaBarrier.csv";
				if (!disorder.load(disorderString, arma::auto_detect) || !ionCoulomb.load(ionCoulombString, arma::auto_detect) 
					|| !DeltaBarrierPotential.load(DeltaBarrierString, arma::auto_detect)) {
					mainout << "-- disorder.csv or ionCoulombAlways.csv not found, aborting ." << endl;
					exit(1);
				}
				else {
					mainout << "total loaded Coulomb ionCoulomb: " << accu(ionCoulomb) << ", total loaded disorder" << accu(disorder) << endl;
					ionCoulomb = ionCoulomb + disorder * disorderStrength + DeltaBarrierPotential;
					disorder.clear();//delete disorder because not used anymore
					DeltaBarrierPotential.clear();
					mainout << "ionCoulomb + disorder * disorderStrength loaded. " << endl;
				}
			}
		}
		else if (externalPotentialGivenFromOutside == true){
			/*something*/
		}
		mainout << "total ionCoulomb: " << accu(ionCoulomb) << endl;
		//calcA_Landau(SeaDensity, BExternalGlobal);
		//calcA_Landau_x(SeaDensity, BExternalGlobal);
		calcA_Landau_x_refined(SeaDensity, BExternalGlobal);
		//mainout << "Landau done: " << accu(A_y_mat) << endl;
		//calcA(SeaDensity, BExternalGlobal); // A_x_mat, A_y_mat are now calculated to be Ax & Ay
		//calcA_Sym(SeaDensity, BExternalGlobal);
		//AddCircumcircleA(useCircumcircleA, BExternalGlobal, whichexternal, ionFilling);
		//writePotentials(); output ranging x in each row, fix y in each row; .csv is not this way
		calcCoulomb(SeaDensity); // VCoulomb or COULOMB::VCoulomb is Now the electron potential instead, DielecX_vs_DielecY already included in green function 1/r definition
		VCoulomb *= scaleHartree;
		mainout << "total VCoulomb Hartree: " << accu(VCoulomb) << endl;		
		// double use of Vcoulomb
		//CalcVxcG(VCoulomb, SeaDensity);
		CalcVxcGradRealFilling(totalPotentialVxc, abs2ToFillingRealB(SeaDensity, BExternalGlobal));
		mainout << "total Vxc: " << accu(totalPotentialVxc) << endl;
		if (continueWithVT) {
			/*this choice has a higher priority*/
			mainout << "-- continueWithVT, oldVT added." << endl;
			mainout << "total oldVT: " << accu(oldVT) << endl;
		}
		else if (continueWithVTfromTable) {
			mainout << "-- continueWithVTfromTable, oldVT added." << endl;
			mainout << "total oldVT: " << accu(oldVT) << endl;
		}
		else { 
			mainout << "-- No oldVT to start with." << endl;
			oldVT.zeros();//otherwise there should be no VT to use, erase in case some weird numbers stored 
		}
		if (verbosity >= Debug && 1 == 2) {
			A_x_mat.save(prefix + "Ax0.csv", arma::csv_ascii);
			A_y_mat.save(prefix + "Ay0.csv", arma::csv_ascii);
			VCoulomb.save(prefix + "VCoulomb0.csv", arma::csv_ascii);
			totalPotentialVxc.save(prefix + "gradVxc0.csv", arma::csv_ascii);
		}
		/* Construct initial hamiltonian.*/
		newSeaDensity = totalPotentialVxc + VCoulomb + ionCoulomb + newVT;//temporarily use newSeaDensity as potential
		double totalPotentialMinimum=FindPotentialMinimum(newSeaDensity,20);//This seems good approximation, strange
		mainout << "totalPotentialMinimum: " << totalPotentialMinimum <<endl;
		//addKinetic_Piers(ham_draft, A_x_mat, A_y_mat);// Mx_vs_My has been considered in kinetic operator function definitions
		addKinetic_Piers_negative(ham_draft, A_x_mat, A_y_mat);
		addPeriodKineticAndOnsite(ham_draft, A_x_mat, A_y_mat, usePer);
		addOnsite(ham_draft, newSeaDensity);
		ham = mapToMMSpmatR<iint, ComplexF64>(Nxy, Nxy, ham_draft); // converted to MMSpmatR struct type, done in SPARSEMAT_HPP in sparsemat.hpps But How? And where is At defined?????
		ham_draft.clear();
		if (verbosity >= Debug) {
			writeHam();// done in main_util.hpp for debug by reading the nonzero elements in ham
			mainout << "--writeHam() Logged zeroth Hamiltonian." << endl;
		}
		auto hamop = create_Hermitian_op<iint>(ham, /*unused*/  Np, /*unused*/ Niter); // done in mklspmat.hpp, for what???? for speed since using MKL?, according to Yang?
		/* Report some initial eneries */
		/*
		 * Compute initial potential energy value, i.e. without V_T*, kinetic term.
		 * TODO: include extrapotential
		 * did yayun copy extrapotential here?
		 */
		mainout << endl << "========================================" << endl << endl;
		/*****************************************
		 * Zeroth iteration.
		 * Run eigensolver with the vanilla conservative subspace size.
		 * Abort as early as possible if eigensolver returns error.
		 * This usually indicates that subspace size (MaxSubspace) and
		 * energy bounds are wildly off.
		 * It is futile to try to search for the correct subspace size automatically.
		 *
		 * Steps:
		 * -- eigensolve
		 * -- calculate energies
		 * -- calculate VT to prime the iteration
		 *
		 * Currently while calibrateEigenBounds can lead to
		 * infinite loop if Eigenlow is too far off.
		 *****************************************/

		iint Nee = 0, lastValidNee, minNee; double estimateDos = 0.0;
		/* This Nsubspace has to do with only nee that can be found; initially use 10 times the particle number as guess; It does not influence geuss subspace TRUE in later iterations */

		/* Diagonalization uses MKL FEAST */
		if (tableNum == LetsStartWhereStops){
			/*yayun 240220, this may cause the zeroth iteration to cost very long time?? */
			mkl_feast_init();// to initialize feast parameters to default values
		}
		mainout << "Zeroth iteration." << endl << endl;
		//addConfinement(ham, confinementHeight); // This should only happen once, so not inside the following while loop: added to ham rather than ham_draft because the latter is removed
		//Nsubspace = min(Nsubspace, MaxSubspace);		
		Nsubspace = min(min(max(max(2000, 8*Np), Nsubspace), MaxSubspace),Nxy);//yayun20240130, add min(?,Nxy)
		if (tableNum==0 && TableAllinOne(tableNum, 10) == -1 /*only purely fresh calculations use this human guess*/) {
			if (knowFeastPara == true) {
				Energylow = knowEnergylow;
				Energyhigh = knowEnergyhigh;
				Nsubspace = knowNsubspace;
				mainout << "use human guess of energy bounds. " << endl;
			}else if (useBoundFromPotential == true) {
				Energylow = totalPotentialMinimum;
				Energyhigh = Energylow + 0.11 * 5 + disorderStrength;//yayun, 231231, before it was "+ 0.11 * 3", now changed to "+ 0.11 * 10" due to small energyhigh when using continueFromfile
				Nsubspace = min(min(max(max(2000, 8*Np), Nsubspace), MaxSubspace),Nxy);//yayun20240130, add min(?,Nxy)
				mainout << "use useBoundFromPotential guess of energy bounds. " << endl;
			}
			lastEmin = totalPotentialMinimum;//otherwise lastEmin is either loaded or calculated; assume lastEmin >= totalPotentialMinimum
		}
		/*
				else if (tableNum>=1 && EfModeFlag){
				fermiEnergy = TableAllinOne(0, 0);
				Energylow = totalPotentialMinimum-0.11*2;
				Energyhigh = fermiEnergy + occuMinimumEkBTFirst * kBT;//this is convenient
				Nsubspace = min(min(1000, 6 * Np), MaxSubspace);
				mainout << "use fixed-Ef guess of energy bounds. " << endl;
			}
		*/
		if (EfModeFlag) {// careful with EfMode and nfix_num issue
			if (tableNum >= 1) { 
				EfMode = 2; 
			fermiEnergy = TableAllinOne(0, 0);//make sure use this fermiEnergy for all the later tableNum
			}
			else { EfMode = 1; }
		}
		//Energylow = -10;Energyhigh = 1;
		mainout <<"Energylow: " << Energylow << ", Energyhigh: " << Energyhigh << endl;
		mainout <<"Nsubspace: " << Nsubspace << ", MaxSubspace: " << MaxSubspace << endl;
		mainout << "a3" << endl;
		int iter = 0; bool checkGo = false;// flag if the diagonalization successfully gives the Ef and ortbials
		do {
			iter = iter + 1;
			mainout << iter << endl;
			if (iter > 10) {
				mainout << "# of searches in zeroth iteration is more than 10. Aborting." << endl;
				exit(1);
			}
			/*compare energy range low to high values*/
			if (Energylow >= Energyhigh) {
				mainout << "You input wrong order of eigenvalue bounds." << endl;
				exit(1);
			}
			/* Diagonlization */
			//assert(Nsubspace <= MaxSubspace); yayun 20240219, this is artificial and not necessary
			if (tableNum == LetsStartWhereStops){
				/*yayun 20240219, to svoid segmentation fault because I assume that with better guess diagonalization requires less storage;  
				That assumption is supported by the fact that resubmit the job goes through zeroth iteration but continuous table for-loop
				does not go through; and the latter takes up more storage due to EigenvecsOld; but later iterations after the Zeroth goes through 
				even with EigenvecsOld defined. The only difference is that later iterations after the Zeroth has good guess == true.
				But useGuessSubspaceZeroth = true here actually does not work since later 2nd zeroth does not go through either.
			    By reducing the storage of EigenvecsOld, later 2nd zeroth go through. 
				What's still not clear is that the time taken for zeroth iteration is as long as before when useGuessSubspaceZeroth = false, which is much 
				longer than later iterations, see E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\80. 
				But this is not always true, see E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\72  */
				useGuessSubspaceZeroth = false;
			}else{
				useGuessSubspaceZeroth = true;
			}			
			TIMEIT_LOG(mainout, Nee = mkl_feast_ev('F', // Full matrix
				ham.n_rows, // Sets the size of the problem %.n_rows	   	number of rows
				(MKL_Complex16*)ham.vals,// Array containing the nonzero elements of either the full matrix A or the upper or lower triangular part of the matrix A, as specified by uplo
				ham.rowstart1,/*Array of length n + 1, containing indices of elements in the array a , such that ia[i] is the index in the array a of the first non-zero element from the row i .The value of the last element ia[n] is equal to the number of nonzeros plus one.*/
				ham.cols1,/*Array containing the column indices for each non-zero element of the matrix A being represented in the array a . Its length is equal to the length of the array a. */
				Energylow, Energyhigh,
				Nsubspace,/* On entry, specifies the initial guess for subspace dimension to be used, 0<m0<=n. Set m0>=m where m is the total number of eigenvalues located in the interval [emin, emax]. If the initial guess is wrong, Extended Eigensolver routines return info=3. */
				/*gives unknown error when Nsubspace is too large than the dimension of the hamiltonian*/
				Eigenvals, reinterpret_cast<MKL_Complex16*>(Eigenvecs.memptr()),
				useGuessSubspaceZeroth)); // Here means FEAST::fpm[4] = false; Do not use guess subspace for zeroth run.		                                      
	//reinterpret_cast convert pointer type without check values
			mainout << "useGuessSubspaceZeroth: " << useGuessSubspaceZeroth << endl;
			if (verbosity >= Debug) {
				reporteigenvalues(-iter, Nee);
				reportOccupation(-iter, Nee);//This value is not assigned in the first run
			}
			/* It's hard to adjust based on a wrong initial guess. So abort. */
			if (FEAST::info != 0) {
				mainout << "FEAST error: " << feast_error_info(FEAST::info) << endl;
				mainout << "Record ees in ees_0.txt." << endl;
				mainout << "Nsubspace : " << Nsubspace << endl;
				reporteigenvaluesFeast(0, Nee);
				mainout << "Aborting." << endl;
				if (FEAST::info == 3){//Warning	Size of the subspace m0 is too small (m0<m)
					mainout << "error 3 with # of eigen state: " << Nee << endl;
					mainout << "error 3 with Eigenvals[Nee-1]: " << Eigenvals[Nee-1] << " Eigenvals[0]: " << Eigenvals[0] << endl;
				}else{
					mainout << "error FEAST::info = " << FEAST::info << ", with # of eigen state: " << Nee << endl;
					mainout << "error FEAST::info = " << FEAST::info << ", with Eigenvals[Nee-1]: " << Eigenvals[Nee-1] << " Eigenvals[0]: " << Eigenvals[0] << endl;
					exit(1);
				}				
			}
			mainout <<"Energylow: " << Energylow << ", Energyhigh: " << Energyhigh << endl;
			mainout << "Previous # of eigen state: " << Nee << "# of Nsubspace: " << Nsubspace << " fixNum: " << fixNum << endl;
			checkGo = false;// ready to modify orbital check
			arma::rowvec sorted_Eigenvals(Nee);
			SortPointerArrayIntoArmaArray(sorted_Eigenvals, Nee, Eigenvals);
			if (calibrateFixEigenSorted(Np, Nee, sorted_Eigenvals, Energyhigh, Energylow, kBT, fermiEnergy, occupationNumber, Nsubspace, estimateDos, fixNum, 
				EfMode, 0/*iter*/, totalPotentialMinimum, lastEmin, nee_wantdecimal) == 0) {
				checkGo = true;
			}
			mainout <<"Energylow: " << Energylow << ", Energyhigh: " << Energyhigh << endl;
			mainout << "# of next Nsubspace: " << Nsubspace  << endl;
			//Nsubspace += 1000;
			Nsubspace = min(max(Nsubspace, int(1.5 * Np) ), MaxSubspace);
			mainout << "# of next Nsubspace: " << Nsubspace  << endl;
		} while (checkGo == false);/*A tedious check defined before main*/
		mainout << "Zeroth iteration passed\n";
		if (EfMode == 1) {
			Nee = fixNum; //number of orbitals used in the following
		}//210210,210326 corrections
		mainout << "a4" << endl;
		setDensityTemperature(newSeaDensity, Eigenvecs, Nee, occupationNumber);
		lastSeaDensity = newSeaDensity; //This is to store previous output eigenstates
		/*
		 * Calculate (instantaneous) DFT functional energy, using immediate eigenvectors.
		 */
		calcEnergyTemperatureRealB(Etotal, KE, Eext, EHartree, Excsum, hamop, newSeaDensity, ionCoulomb, Eigenvecs, Nee, occupationNumber, BExternalGlobal);
		double oldKE = KE, oldNp = Np;
		mainout << "Total Energy: " << Etotal << " = " << KE << "(KE)+" << Eext << "(Eext)+" << EHartree << "(EHartree)+" << Excsum << "(Exc)." << endl;
		recordEnergy(energyFile + TableRoundNum, 0, Etotal, KE, Eext, EHartree, Excsum);

		if (!continueWithVT && !continueWithVTfromTable)//otherwise stick to whatever is loaded to use as initial VT, which is modified only in later iterations(iter>=2)
		{
			// calculate VT to prime the iteration only if none to start with
			calcDynMomentumTemperature(dynP_x, dynP_y, newSeaDensity, Eigenvecs, Nee, occupationNumber, A_x_mat, A_y_mat);//A_x_mat, A_y_mat are not updated, only dynP_x, dynP_y updated
			calcVT(dynP_x, dynP_y);// using namespace VTSTAR; makes VT in scope here
			//oldVT = (VT_x_mat + VT_y_mat) * scaleVT *0;/*this is very important, means turn on VT slowly; Yayun0402*/
			oldVT = (VT_x_mat + VT_y_mat) * scaleVT;/*turn on VT immediately; Yayun0410*/
			oldVT.save(prefix + "VT0.csv", arma::csv_ascii); //write this freshingly calculated VT only when no VT0con.csv. This means switch on VT suddenly
			//also adjust energy bounds shift by VT, because switch on VT suddenly
			Energyhigh += oldVT.max();
			Energylow += oldVT.min();
			Nsubspace += estimateDos * (oldVT.max() - oldVT.min());
			lastVT = oldVT;
		}

		if (verbosity >= Info) {
			reporteigenvalues(0, Nee);
			reportOccupation(0, Nee);//This value is now assigned
		}
		//{	//to change MaxSubspace
			//this brackets is very important because temperatory variables defined inside will be released outside
			lastValidNee = Nee;
			if (tableNum == LetsStartWhereStops){
				{
				mainout << "change MaxSubspace" << endl;
				MaxSubspace = min(max(max(10 * Np, 2 * Nsubspace),2000), MaxSubspace); 
				//Yayun 240219. to avoid segmentation fault due to machine storage limit
				Nee = min(Nee, MaxSubspace);//Yayun 240219. to avoid segmentation fault during copy, logically necessary
				//define and copy data to Old arrays after size determined
				occupationNumberOld = new double[MaxSubspace];
				EigenvalsOld = new double[MaxSubspace];
				std::copy(Eigenvals, Eigenvals + Nee, EigenvalsOld);
				std::copy(occupationNumber, occupationNumber + Nee, occupationNumberOld);
				//update working arrays to new size and refill its data
				Eigenvecs.resize(Nxy, MaxSubspace);
				EigenvecsOld = Eigenvecs; //used to determine updateratio through overlap in function "evsOverlapDiff"
				EigenvecsOld.resize(Nxy, Nee);
				Eigenvals = new double[MaxSubspace];
				occupationNumber = new double[MaxSubspace];
				//refill
				std::copy(EigenvalsOld, EigenvalsOld + Nee, Eigenvals);
				std::copy(occupationNumberOld, occupationNumberOld + Nee, occupationNumber);
				//delete[] Eigenvals; 
				Nsubspace = min(Nsubspace, MaxSubspace);
				mainout << "Maximum eigensubspace size, MaxSubspace, is reset as " << MaxSubspace << endl << endl;
				mainout << "Nsubspace is: " << Nsubspace << endl << endl;
				mainout << "# of eigen state is reduced to or remains the same as: " << Nee << endl;
				}
			}
		//}
		recordoccuAllbinary(occuAllFilebin + TableRoundNum, 0, Nee);
		recordEigenAllbinary(eigenAllFilebin + TableRoundNum, 0, Nee);

		double occudiff = 0.0, eesdiff = 0.0, vectorDiff = 0.0;
		arma::rowvec evsdiff(MaxSubspace, arma::fill::zeros);
		arma::cx_rowvec evsoverlap(MaxSubspace, arma::fill::zeros);
		arma::urowvec connection(1, MaxSubspace);
		//for mixing rate manipulation
		arma::rowvec occudeltaPre(MaxSubspace, arma::fill::zeros), occudelta(MaxSubspace, arma::fill::zeros), eesdeltaPre(MaxSubspace, arma::fill::zeros), eesdelta(MaxSubspace, arma::fill::zeros);
		double rateOccu = 0.0, thetaOccu = 0.0, expansionOccu = 0.0, deltaOccu = 0.0, eesrateE = 0.0, eesthetaE = 0.0, eesexpansionE = 0.0, eesdeltaE = 0.0;// for mixing rate manipulation
		int minNeePre = 0, Bigiter = 0;
		double BigPeriodDiff = 0.0;// check big period;
		arma::mat SeaDensityBig;
		bool ztoolargelast = false;
		/*to start with a new table parameter, prepare update ratios to conservative values ~ 0*/
		//UpdateRatio = UpdateRatioLowBound;UpdateRatioVT = UpdateRatio; yayun20240219

		recordfermi(fermiFile+ TableRoundNum, 0, lastValidNee, lastValidNee, Nsubspace, Energyhigh, Energylow, fermiEnergy, occudiff, eesdiff, relVTDiff, relVTDiffAim, relDensityDiff, relDensityDiffAim, vectorDiff,
			UpdateRatio, rateOccu, thetaOccu, expansionOccu, deltaOccu, BigPeriodDiff, eesrateE, eesthetaE, eesexpansionE, eesdeltaE, totalPotentialMinimum);

		/*****************************************
		 * iteration
		 * Steps:
		 * -- Update electron density.
		 * -- Update Hamiltonian / energy functional.
		 * -- Eigensolve, while updating energy bounds
		 * ---- to ensure right number of eigenstates.
		 * -- Calculate energies.
		 * -- Reduce subspace size if possible.
		 * -- Test stoping criteria, ie. density change threshold.
		 *****************************************/

		mainout << endl << "========================================" << endl << endl;
		mainout << "DFT iteration." << endl << endl;


		int save_cycle1 = 0, save_cycle2 = 0;//initialize save_cycle, which is used to cycle between files every save_cycleNum period, this is done in case storage corrupt
        //above should be initialized here for each new round of iteration

		for (auto iter = 1; iter <= Niter; iter++) {
			mainout << endl << "----------------------------------------" << endl;
			mainout << "Recursion step: " << iter << endl;
			mainout << "Nsubspace: " << Nsubspace <<" Np: " << Np << endl;
			//assert(Nsubspace >= Np); yayun 20240219, this is artificial and not necessary
			 //Update electron density
			 //UpdateRatio = 1.0;
			SeaDensity = newSeaDensity * UpdateRatio + SeaDensity * (1 - UpdateRatio);
			/* Update Hamiltonian / energy functional.*/			 
			//calcA_Landau(SeaDensity, BExternalGlobal);
			//calcA_Landau_x(SeaDensity, BExternalGlobal);
			mainout << "accu(SeaDensity): " << accu(SeaDensity) << endl;
			calcA_Landau_x_refined(SeaDensity, BExternalGlobal);
			//calcA_Sym(SeaDensity, BExternalGlobal);
			//calcA(SeaDensity, BExternalGlobal); // A_x_mat, A_y_mat are now calculated to be Ax & Ay
			//AddCircumcircleA(useCircumcircleA, BExternalGlobal, whichexternal, ionFilling);
			calcCoulomb(SeaDensity); // VCoulomb or COULOMB::VCoulomb is Now the electron potential instead, DielecX_vs_DielecY already included in green function 1/r definition
			VCoulomb *= scaleHartree;
			//mainout << "accu(VCoulomb): " << accu(VCoulomb) << endl;
			// calculate V_T^*
			calcDynMomentumTemperature(dynP_x, dynP_y, SeaDensity, Eigenvecs, Nee, occupationNumber, A_x_mat, A_y_mat);
			calcVT(dynP_x, dynP_y);
			newVT = (VT_x_mat + VT_y_mat) * scaleVT;
			//mainout << "accu(newVT): " << accu(newVT) << endl;
			{
				relVTDiff = accu(abs(newVT - oldVT)) / accu(abs(newVT - newVT.min()));//mainout << "VT change summed over all space: " << relVTDiff << endl;
				relVTDiffAim = accu(abs(newVT - lastVT)) / accu(abs(newVT - newVT.min()));
			}
			lastVT = newVT;//fresh output VT
			newVT = newVT * UpdateRatioVT + oldVT * (1 - UpdateRatioVT);
			/*{
				INFO(newVT.max());INFO(newVT.min());
				calcVT(XX,YY);
				auto p = (VT_x_mat + VT_y_mat);
				INFO(p.max());INFO(p.min());
				//	assert(p.max() < 1e-6);
				//	assert(p.min() >-1e-6);
			}*/
			//mainout << "-- A & Vcoulomb convolved." << endl;
			oldVT = newVT;//VT used in the coming diagonalization 
			CalcVxcGradRealFilling(totalPotentialVxc, abs2ToFillingRealB(SeaDensity, BExternalGlobal));
			// rebuild ham
// principles when adding terms: values of variables intact, rescale considered by multiplying coefficients
			newSeaDensity = totalPotentialVxc + VCoulomb + ionCoulomb + newVT;//temporarily use newSeaDensity as potential
			//mainout << "accu(newSeaDensity): " << accu(newSeaDensity) << endl;
			totalPotentialMinimum = FindPotentialMinimum(newSeaDensity, 20);
			std::fill(ham.vals, ham.vals + ham.nnz, (ComplexF64)(0));//initialize the sparse matrix ham
			//addKinetic_Piers(ham, A_x_mat, A_y_mat);
			//mainout << "accu(A_x_mat): " << accu(A_x_mat) << endl;
			//mainout << "accu(A_y_mat): " << accu(A_y_mat) << endl;
			addKinetic_Piers_negative(ham, A_x_mat, A_y_mat);
			addPeriodKineticAndOnsite(ham, A_x_mat, A_y_mat, usePer);//not used
			//mainout << "totalPotentialMinimum: " << totalPotentialMinimum << endl;
			addOnsite(ham, newSeaDensity);
			//addConfinement(ham, confinementHeight);
			mainout << "-- Hamiltonian assembled." << endl;

			//there used to be a weird error that eigenvalues are not found in first iteration. but if the following 3 lines are moved to later then fine? Why??? Yayun    
			//calcEnergyTemperature(Etotal, KE, Eext, EHartree, Excsum, hamop, newSeaDensity, ionCoulomb, Eigenvecs, Nee, occupationNumber);
			//mainout << "Total Energy: " << Etotal << " = " << KE << "(KE)+" << Eext << "(Eext)+" << EHartree << "(EHartree)+" << Excsum << "(Exc)." << endl;
			/*
			 * Eigensolve, while updating energy bounds.
			 * If subspace size too small, we double it.
			 */
			mainout << "a3" << endl;
			do {
				mainout << "Energylow: " << Energylow << ", Energyhigh: " << Energyhigh << endl;
				/* Diagonalization */
				//assert(Nsubspace <= MaxSubspace); yayun 20240219, this is artificial and not necessary
				TIMEIT_LOG(mainout, Nee = mkl_feast_ev('F', 
				ham.n_rows, (MKL_Complex16*)ham.vals, ham.rowstart1, ham.cols1, 
				Energylow, Energyhigh, Nsubspace, Eigenvals, 
				reinterpret_cast<MKL_Complex16*>(Eigenvecs.memptr()), true));
				mainout << "Energylow: " << Energylow << ", Energyhigh: " << Energyhigh << endl;
				if (FEAST::info != 0) {
					mainout << "FEAST error: " << feast_error_info(FEAST::info) << endl;
					if (FEAST::info == 3) {
						mainout << "Nsubspace too small: " << Nsubspace << endl;
					}
					else {
						mainout << "FEAST error != 3 occur. " << endl;
						reporteigenvaluesFeast(iter, Nee);
						Nee = -1;
					}
				}
				//INFO_LOG(mainout,Nee);
				//INFO_LOG(mainout, Nsubspace);//# of next Nsubspace
				//INFO_LOG(mainout, fixNum);
				mainout << "Previous Nee: " << Nee << ",Nsubspace: " << Nsubspace << ",fixNum: " << fixNum << endl;
				checkGo = false;// ready to modify orbital check
				arma::rowvec sorted_Eigenvals(Nee);
				SortPointerArrayIntoArmaArray(sorted_Eigenvals, Nee, Eigenvals);
				if (calibrateFixEigenSorted(Np, Nee, sorted_Eigenvals, Energyhigh, Energylow, kBT, fermiEnergy, occupationNumber, Nsubspace, estimateDos, fixNum, 
					EfMode, iter, totalPotentialMinimum, lastEmin, nee_wantdecimal) == 0) {
					/* calibrateFixEigenZeroth(Np, Nee, Eigenvals, Energyhigh, Energylow, kBT, fermiEnergy, occupationNumber, Nsubspace, estimateDos, fixNum, 
					EfMode, iter, totalPotentialMinimum, lastEmin, nee_wantdecimal) == 0 */
					checkGo = true;
				}
				mainout << "# of next Nsubspace: " << Nsubspace  << endl;
				//mainout << "Current Nee: " << Nee << ",Nsubspace: " << Nsubspace << ",fixNum: " << fixNum << endl;
				Nsubspace = min(max(Nsubspace, int(1.5 * Np) ), MaxSubspace);
				/*yayun, 20240127, it happens that Nsubspace is too small when FERMI_ENERGY_STATUS_NUM == 4 "E:\sharedata01\yhu\workspace\data\CFDFTsquare\LXfracFixNee240116\28\slurm-35187_1.out"	
				but this should not have been a problem or I should say that it is an artificial problem. So at least fix by hand.*/
				mainout << "# of next Nsubspace: " << Nsubspace  << endl;
			} while (checkGo == false);
			mainout << "a4" << endl;
			TableAllinOne(tableNum, 22) = Nee;//22number of orbitals really get
			if (EfMode == 1) {
				Nee = fixNum;
			}//210210,210326 corrections
			//} while( calibrateEigenBoundsTemperature(Np, Nee, Eigenvals, Energyhigh, Energylow, kBT, fermiEnergy, occupationNumber, Nsubspace, estimateDos) != 0 );
			//setDensity(newSeaDensity, Eigenvecs, Np);
			setDensityTemperature(newSeaDensity, Eigenvecs, Nee, occupationNumber);
			/* DFT functional energy calculation */
			//calcEnergyTemperature(Etotal, KE, Eext, EHartree, Excsum, hamop, newSeaDensity, ionCoulomb, Eigenvecs, Nee, occupationNumber);
			calcEnergyTemperatureRealB(Etotal, KE, Eext, EHartree, Excsum, hamop, newSeaDensity, ionCoulomb, Eigenvecs, Nee, occupationNumber, BExternalGlobal);
			//mainout << "Total Energy: " << Etotal << " = " << KE << "(KE)+" << Eext << "(Eext)+" << EHartree << "(EHartree)+" << Excsum << "(Exc)." << endl;
			//record total energy, occu and eigen
			recordEnergy(energyFile + TableRoundNum, iter, Etotal, KE, Eext, EHartree, Excsum);
			recordoccuAllbinary(occuAllFilebin + TableRoundNum, iter, Nee);
			//recordEigenAll(eigenAllFile, iter, Nee);	
			recordEigenAllbinary(eigenAllFilebin + TableRoundNum, iter, Nee);
			//compare with previous  eig, occu  and eigenvector
			minNee = min(lastValidNee, Nee);
			occupationDiff(minNee, Nee, occudiff, occudelta);
			eesDiff(minNee, Nee, eesdiff, eesdelta);//old ees and occu are also updated here.
			evsOverlapDiff(minNee, Nee, evsdiff, evsoverlap, vectorDiff);//particle-equivalennce change in each single-particle orbital
			//save into files
			recordEigendiffbinary(evsdiffFilebin + TableRoundNum, iter, minNee, evsdiff);
			recordEigenoverlapbinary(evsoverlapFilebin + TableRoundNum, iter, minNee, evsoverlap);
			
			//mixing rate manipulation
			if (iter == 1) {
				occudeltaPre.cols(0, minNee - 1) = occudelta.cols(0, minNee - 1);
				eesdeltaPre.cols(0, minNee - 1) = eesdelta.cols(0, minNee - 1);
				minNeePre = minNee;
				ztoolargelast = false;
			}
			else {
				arma::rowvec A_Vec_InPlace(occupationNumber, Nee, false, true);//convert occupationNumber to arma
				newUpdateRate(ztoolargelast, minNeePre, minNee, UpdateRatioTemp, rateOccu, thetaOccu, expansionOccu, deltaOccu/*the ratio of change over occudeltaPre*/,
					occudeltaPre, occudelta, evsdiff, A_Vec_InPlace, relVTDiff, eesdeltaPre, eesdelta, eesrateE, eesthetaE, eesexpansionE, eesdeltaE);
				//yayun20240219: now UpdateRatio is set in syspara.hpp	, here UpdateRatio will be used along the table once set
				if (iter>30){
					UpdateRatio = 0.02;	
				}else if(iter>80){
					UpdateRatio = 0.01;	
				}else if(iter>150){
					UpdateRatio = 0.005;	
				}else if(iter>300){
					UpdateRatio = 0.001;	
				}
				//UpdateRatio = 0.05;//yayun240130							
				UpdateRatioVT = UpdateRatio;
			}
			/*
			if(verbosity >= Info && iter % eesfrequency == 1) {
				reporteigenvalues(iter, Nee);
				reportOccupation(iter, Nee);
			}
			if (verbosity >= Info && iter % evsfrequency == 1) {
				reporteigenvec(iter, Nee);
			}*/
			relDensityDiff = accu(abs(newSeaDensity - SeaDensity)) / Np;
			relDensityDiffAim = accu(abs(newSeaDensity - lastSeaDensity)) / Np;
			//mainout << "Density/Np change summed over all space: "<<relDensityDiff<<endl;
			mainout << "relDensityDiff: " << relDensityDiff << endl;
			//big period tracking
			if (relDensityDiff < enterBigThreshold && relDensityDiffAim < enterBigThreshold && relVTDiff < enterBigThreshold 
				&& relVTDiffAim < enterBigThreshold && Bigiter == 0) {
				Bigiter = iter;
				BigPeriodDiff = Bigiter;//the first non-zero element of BigPeriodDiff is the reference density for comparison
				SeaDensityBig = SeaDensity;
				TableAllinOne(tableNum, 13) = Bigiter;//13RhoBigNum
			}
			if (Bigiter > 0) {
				BigPeriodDiff = accu(abs(SeaDensityBig - SeaDensity));//effective from step Bigiter+1, because Bigiter defines SeaDensityBig
			}
			/*update data convergence considerations*/
			AllbigDiff(iter % lookBackStepNum) = BigPeriodDiff;
			AlldensityDiff(iter % lookBackStepNum) = relDensityDiff;
			AllVTDiff(iter % lookBackStepNum) = relVTDiff;
			AllkineticEnergy(iter % lookBackStepNum) = KE - oldKE;
			Allnewparticle(iter % lookBackStepNum) = Np-oldNp;
			oldKE = KE; oldNp = Np;

			recordfermi(fermiFile + TableRoundNum, iter, Nee, minNee, Nsubspace, Energyhigh, Energylow, fermiEnergy, occudiff, eesdiff, 
				relVTDiff, relVTDiffAim, relDensityDiff, relDensityDiffAim, vectorDiff,
				UpdateRatio, rateOccu, thetaOccu, expansionOccu, deltaOccu, BigPeriodDiff, eesrateE, eesthetaE, eesexpansionE, eesdeltaE,
				totalPotentialMinimum);
			lastValidNee = Nee;
			if (verbosity >= Debug && ((iter % denfrequency == 1) || (iter == 5 || iter == 10 || iter == 50))) {
				/*save at the same frequency as density, update only when saving*/
				saveHamToFile(prefix, TableRoundNum, save_cycle1);
				saveDataToFile(prefix, TableRoundNum, save_cycle1, SeaDensity, newSeaDensity, oldVT, A_x_mat, A_y_mat, VCoulomb, totalPotentialVxc);
				assignTableAllinOne(TableAllinOne, tableNum, fermiEnergy/*1*/, Etotal/*2*/, KE/*3*/, EHartree/*4*/, Eext/*5*/, Excsum/*6*/,
					iter/*7*/, relDensityDiff/*8*/, relDensityDiffAim/*9*/, save_cycle1/*10*/, Energylow/*11*/, Energyhigh/*12*/, fixNum/*13*/, Nsubspace/*14*/,
					totalPotentialMinimum/*15*/, lastEmin/*16*/, nee_wantdecimal/*17*/, Nee/*18*/);
			    TableAllinOne.save(prefix + "TableAllinOne.csv", arma::csv_ascii);			
			    TableAllinOne.save(prefix + "TableAllinOne" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			    /*
			    AllbigDiff.save(prefix + "AllbigDiff" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);;
			    AlldensityDiff(iter % lookBackStepNum) = relDensityDiff;
			    AllVTDiff(iter % lookBackStepNum) = relVTDiff;
			    AllkineticEnergy(iter % lookBackStepNum) = KE - oldKE;
			    Allnewparticle(iter % lookBackStepNum) = Np - oldNp;
			    */
            save_cycle1 += 1;
            save_cycle1 = save_cycle1 % save_cycleNum1;
			}
			//bool whataboutVT = AllVTDiff.max() < relDensityThreshold;
			//mainout << "whataboutVT1 : " <<  whataboutVT  << endl;
			//whataboutVT = abs(AllkineticEnergy.max() - AllkineticEnergy.min()) * 0.0 < abs(KineticThreshold * KE);
			//mainout << "whataboutVT2 : " << whataboutVT << endl;
			//mainout << "AllVTDiff.max() : " << AllVTDiff.max()  << endl;
			//mainout << "AllkineticEnergy.min() AllkineticEnergy.max() : " << AllkineticEnergy.min() << "  " << AllkineticEnergy.max() << endl;
			if (relDensityDiff < relDensityThreshold && ConvergenceRecentNoBig(iter, Bigiter, KE, AllbigDiff, AlldensityDiff, AllVTDiff, AllkineticEnergy, Allnewparticle)) {//
				saveHamToFile(prefix, TableRoundNum, save_cycle1);
				saveDataToFile(prefix, TableRoundNum, save_cycle1, SeaDensity, newSeaDensity, oldVT, A_x_mat, A_y_mat, VCoulomb, totalPotentialVxc);
				assignTableAllinOne(TableAllinOne, tableNum, fermiEnergy/*1*/, Etotal/*2*/, KE/*3*/, EHartree/*4*/, Eext/*5*/, Excsum/*6*/,
					iter/*7*/, relDensityDiff/*8*/, relDensityDiffAim/*9*/, save_cycle1/*10*/, Energylow/*11*/, Energyhigh/*12*/, fixNum/*13*/, Nsubspace/*14*/,
					totalPotentialMinimum/*15*/, lastEmin/*16*/, nee_wantdecimal/*17*/, Nee/*18*/);				
				recordoccuConvergeOrFinished(convergedOccu, tableNum, Nee);
				recordEigenConvergeOrFinished(convergedEigen, tableNum, Nee);
				reporteigenvecFinished(prefix + "evs" + TableRoundNum + to_string((int)save_cycle1) + ".bin", Nee);
				/*save above because this iter may not have been saved*/
				TableAllinOne(TemperatureNum * MySecondSetNum, 0)=tableNum;//0tableNumFinished
				TableAllinOne(tableNum, 18) = iter;//18iterConverge, it also flags convergence
				TableAllinOne(tableNum, 19) = 1;//19ConvergeFlag, it is purposely used to flag convergence; Converged=1; has been started=0; finished but not converged=-2;default=-1;
				TableAllinOne.save(prefix + "TableAllinOne.csv", arma::csv_ascii);
				TableAllinOne.save(prefix + "TableAllinOne" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
                //TableAllinOne = TableAllinOne * (-1);//initialize to -1;
/*0FermiEnergy,1TotalAngualr,2TotalEigenEnergy,3TotalKinetic,4TotalKohnShamPotentialEnergy,5TotalDensityEnergy,6ElectronHartreeEnergy,7DiskPotentialEnergy,
8KineticDensityEnergy,9Exchange-correlation-Energy,10iterationNum,11OverlapRatio,12OverlapRatioNum,13RhoBigNum,14roundNum,15savecycle1,16savecycle2,17tableNumStartStage,
18iterConverge,19ConvergeFlag,20Energylow,21Energyhigh*/
/*TableAllinOne(TemperatureNum* MySecondSetNum, );0tableNumFinished,1tableNumStarted,2notConverged,3savedB,4eratureNum,5MySecondSetNum,
6save_cycleNum1 for density,7save_cycleNum2 for potentials*/ 
				mainout << "Density change smaller than threshold: " << relDensityThreshold << '.' << endl
					<< "Finishing." << endl;
				break;
			}
			lastSeaDensity = newSeaDensity;//lastSeaDensity is the previous newSeaDensity
			/*
			if (iter % 10 == 3) {
				int continueORnot;
				mainout << "continue1 or not2: " << endl;
				std::cin >> continueORnot;
				if (continueORnot == 1) {
					mainout << "continued." << endl;
				}
				else if(continueORnot == 2) {
					mainout << "not continued." << endl;
					exit(-1);
				}
			}
			*/
		}//closure of "for (auto iter = 1; iter <= Niter; iter++) {"
		if (TableAllinOne(tableNum, 19) != 1) {//when table is finished but not converged
			saveHamToFile(prefix, TableRoundNum, save_cycle1);
			saveDataToFile(prefix, TableRoundNum, save_cycle1, SeaDensity, newSeaDensity, oldVT, A_x_mat, A_y_mat, VCoulomb, totalPotentialVxc);
			assignTableAllinOne(TableAllinOne, tableNum, fermiEnergy/*1*/, Etotal/*2*/, KE/*3*/, EHartree/*4*/, Eext/*5*/, Excsum/*6*/,
				iter/*7*/, relDensityDiff/*8*/, relDensityDiffAim/*9*/, save_cycle1/*10*/, Energylow/*11*/, Energyhigh/*12*/, fixNum/*13*/, Nsubspace/*14*/,
				totalPotentialMinimum/*15*/, lastEmin/*16*/, nee_wantdecimal/*17*/, Nee/*18*/);
			recordoccuConvergeOrFinished(convergedOccu, tableNum, Nee);
			recordEigenConvergeOrFinished(convergedEigen, tableNum, Nee);
			reporteigenvecFinished(prefix + "evs" + TableRoundNum + to_string((int)save_cycle1) + ".bin", Nee);
			/*save above because this (which is the last) iter may not have been saved*/
			TableAllinOne(TemperatureNum* MySecondSetNum, 0) = tableNum;//0tableNumFinished
			TableAllinOne(tableNum, 19) = -2;//19ConvergeFlag, it is purposely used to flag convergence; Converged=1; has been started=0; finished but not converged=-2;default=-1;
			TableAllinOne.save(prefix + "TableAllinOne.csv", arma::csv_ascii);
			TableAllinOne.save(prefix + "TableAllinOne" + TableRoundNum + to_string((int)save_cycle1) + ".csv", arma::csv_ascii);
			mainout << "Density not converged, but iteration finishes here." << endl;
		}

	}
		mainout << endl << "========================================" << endl;
		mainout << "DFT ended" << endl;
		{
			auto t1 = high_resolution_clock::now();
			const double dt = duration_cast<seconds>(t1 - startTime).count();
			mainout << "Total elapsed time " << dt << "s." << endl;
		}
		mainout << "========================================" << endl;
		return 0;
}

