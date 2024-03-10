// This file is under public domain.
#ifndef ELECTROMAG_HPP
#define ELECTROMAG_HPP

/*****************************************
 * Classical electromagnetic potential generators,
 * and V_T^* term that is the variation of A of wave functions
 * using MKL vsl convolution utilities.
 *****************************************/

#include "commoninclude.hpp"
#include "syspara.hpp"
#include <armadillo>
#include <mkl_types.h>
#include <mkl_vsl.h>

/*
 * mkl convolution
 * They use a "open boundary" type convention which only allow
 * summands t>=0 to be nonzero. In other words, all beyond the
 * boundaries are cutoff.
 */
/* <mkl.h> */
/* https://software.intel.com/en-us/mkl-developer-reference-c-convolution-and-correlation */
// auto status = int vslzConvNewTaskX(task, mode, dims, xshape, yshape, zshape, x, xstride);
// auto status = int vslsConvExecX(VSLConvTaskPtr task, const MKL_Complex16[] y, const int[] ystride, const MKL_Complex16[] z, const int[] zstride);

void vslerror(const MKL_INT status, const char *msg) {
	using namespace std;
	if( status != VSL_STATUS_OK) {
		cout << msg << ": " << status << " of " << VSL_STATUS_OK << endl;
		exit(status);
	}
}

/*
 * centerGreen2D
 * Self-contained convolution task structure. Not used in the program yet.
 * assuming continuous storage, since arma
 */
typedef struct centerGreen2D {
	MKL_INT yzshape[2];
	arma::mat Green;
	VSLConvTaskPtr task;
	void createTask() { /* there is no delete */
		auto mode = VSL_CONV_MODE_AUTO;
		auto gf = reinterpret_cast<const double*>(Green.memptr());
		MKL_INT dims = 2;
		MKL_INT xshape[2] = {(MKL_INT)Green.n_rows,(MKL_INT)Green.n_cols};
		MKL_INT xstride[2] = {1,xshape[0]};
		MKL_INT start[2] = {(MKL_INT)(Green.n_rows/2),(MKL_INT)(Green.n_cols/2)};

		auto status = vsldConvNewTaskX(&task, mode, dims, xshape, yzshape, yzshape, gf, xstride);
		vslerror(status, "Failed vsldConvNewTaskX");
		status = vslConvSetStart(task, start);
		vslerror(status, "Failed vsldConvSetStart");
	}
	void convIt(const arma::mat &input, arma::mat &output) {
		auto vc = reinterpret_cast<double*>(output.memptr());
		auto rho = reinterpret_cast<const double*>(input.memptr());
		MKL_INT ystride[2] = {1,(MKL_INT)input.n_rows};

		if(input.n_cols != yzshape[1] || input.n_rows != yzshape[0]) {
			std::cerr << "Wrong input matrix shape for convolution." <<std::endl;
			std::exit(1);
		}
		if(output.n_cols != yzshape[1] || output.n_rows != yzshape[0]) {
			std::cerr << "Wrong output matrix shape for convolution." <<std::endl;
			std::exit(1);
		}
		
		auto status = vsldConvExecX(task,rho,ystride,vc,ystride);
		vslerror(status, "Failed vsldConvExecX");
	}
} centerGreen2D;

/* Good example initializer for centerGreen2D. Not used as much. */
template<typename iint, typename T>
centerGreen2D& buildCenterGreen2D(iint nx, iint ny, T initFn) {
	centerGreen2D *g = new centerGreen2D;
	initFn(g->Green);
	g->yzshape[0] = nx;
	g->yzshape[0] = ny;
	g->createTask();
	return *g;
}

/*****************************************
 * Vector potential
 *****************************************/
namespace VECTORPOT {//vector potential
	arma::mat greenA_x, greenA_y; // Green's functions.
	arma::mat A_x_mat, A_y_mat; // Ax & Ay. Convolution results.
	/* MKL_VSL */
	VSLConvTaskPtr taskA_x, taskA_y;
};


void initGreenA() {
	using namespace VECTORPOT;

	/* filling Green function matrices */
	auto ayl = 2*Ny-1;
	auto axl = 2*Nx-1;
	auto mx=(axl-1)/2, my=(ayl-1)/2;
	greenA_x = arma::mat(axl,ayl);
	greenA_y = arma::mat(axl,ayl);
	
	const auto c = ax*ay/2/Pi;

	for(iint j=0; j<ayl; j++) {
		for(iint i=0; i<axl; i++) {
			double r2 = (i-mx)*(i-mx)*ax*ax + (j-my)*(j-my)*ay*ay;//in units of (l_B)^2
			greenA_x(i,j)=c*(my-j)*ay/r2;
			greenA_y(i,j)=c*(i-mx)*ax/r2;
		}
	}
	greenA_x(mx,my) = 0;
	greenA_y(mx,my) = 0;
	
	A_x_mat = arma::mat(Nx,Ny);
	A_y_mat = arma::mat(Nx,Ny);
}

void createConvTaskA() {
	using namespace VECTORPOT;

	auto mode = VSL_CONV_MODE_AUTO;
	auto gA_x = reinterpret_cast<const double*>(greenA_x.memptr());
	auto gA_y = reinterpret_cast<const double*>(greenA_y.memptr());
	MKL_INT dims = 2;
	MKL_INT xshape[2] = {(MKL_INT)greenA_x.n_rows,(MKL_INT)greenA_x.n_cols};
	MKL_INT yshape[2] = {(MKL_INT)A_x_mat.n_rows,(MKL_INT)A_x_mat.n_cols};
	MKL_INT xstride[2] = {1,xshape[0]};
    /* The default start is {0,0,...}, 0-based. */
	MKL_INT start[2] = {(MKL_INT)(greenA_x.n_rows/2),(MKL_INT)(greenA_x.n_cols/2)};

	auto status = vsldConvNewTaskX(&taskA_x, mode, dims, xshape, yshape, yshape, gA_x, xstride);
	vslerror(status, "Failed vsldConvNewTaskX for A_x");
	status = vslConvSetStart(taskA_x, start);
	vslerror(status, "Failed vslConvSetStart for A_x");
	
	status = vsldConvNewTaskX(&taskA_y, mode, dims, xshape, yshape, yshape, gA_y, xstride);
	vslerror(status, "Failed vsldConvNewTaskX for A_y");
	status = vslConvSetStart(taskA_y, start);
	vslerror(status, "Failed vslConvSetStart for A_y");
}

/* Convolve to get A */
void calcA(const arma::mat &density /*wave function squared*/, const double& BExternalGlobaltemp) {
	using namespace VECTORPOT;

	MKL_INT ystride[2] = {1,static_cast<int>(density.n_rows)};// Why define a MKL_INT?

	arma::mat B_reduced = (1.0 + BExternalGlobaltemp) - 2.0 * abs2ToFilling(density);
	
	auto rho = reinterpret_cast<const double*>(B_reduced.memptr());
	auto A_x = reinterpret_cast<double*>(A_x_mat.memptr());
	auto A_y = reinterpret_cast<double*>(A_y_mat.memptr());

	auto status = vsldConvExecX(taskA_x,rho,ystride,A_x,ystride);
	vslerror(status, "Failed vsldConvExecX for A_x");

	status = vsldConvExecX(taskA_y,rho,ystride,A_y,ystride);
	vslerror(status, "Failed vsldConvExecX for A_y");
}

/* Convolve to get A */
void calcA_Sym(const arma::mat &density /*wave function squared*/, const double& BExternalGlobaltemp) {
	using namespace VECTORPOT;

	MKL_INT ystride[2] = {1,static_cast<int>(density.n_rows)};// Why define a MKL_INT?

	arma::mat B_reduced = (1.0 + BExternalGlobaltemp) -0.0*abs2ToFilling(density);
	
	auto rho = reinterpret_cast<const double*>(B_reduced.memptr());
	auto A_x = reinterpret_cast<double*>(A_x_mat.memptr());
	auto A_y = reinterpret_cast<double*>(A_y_mat.memptr());

	auto status = vsldConvExecX(taskA_x,rho,ystride,A_x,ystride);
	vslerror(status, "Failed vsldConvExecX for A_x");

	status = vsldConvExecX(taskA_y,rho,ystride,A_y,ystride);
	vslerror(status, "Failed vsldConvExecX for A_y");
}


void calcA_Landau(const arma::mat &density /*wave function squared*/, const double& BExternalGlobaltemp) {
	using namespace VECTORPOT;//yayun, 240101, sure necessary to reuse but not sure how to avoid reuse
	
	//arma::mat B_reduced = (1.0 + BExternalGlobaltemp) -2.0*abs2ToFilling(density);
	arma::mat B = (1.0 + BExternalGlobaltemp) -0.0*abs2ToFilling(density);
	
	arma::mat A_x=0.0*abs2ToFilling(density);
	A_x.zeros();
	arma::mat A_y=A_x;	
	
	auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
	auto x_origin = round(mx);
	for (iint j = 0; j < Ny; j++) {
		A_y(x_origin, j) = 0.0;
	}
	
    for (iint i = x_origin + 1; i < Nx; i++) {
		for (iint j = 0; j < Ny; j++) {
			A_y(i, j) = A_y(i-1, j) + B(i-1, j) * ax;
		}
	}
	
    for (iint i = x_origin - 1; i >= 0; i--) {
		for (iint j = 0; j < Ny; j++) {
			A_y(i, j) = A_y(i+1, j) - B(i, j) * ax;
		}
	}
	
	//cout << "debug: " << accu(A_x) << endl;

	A_x_mat = A_x;
	
	A_y_mat = A_y;
	
}

void calcA_Landau_x(const arma::mat &density /*wave function squared*/, const double& BExternalGlobaltemp) {
	using namespace VECTORPOT;//yayun, 240101, sure necessary to reuse but not sure how to avoid reuse
	
	//arma::mat B_reduced = (1.0 + BExternalGlobaltemp) -2.0*abs2ToFilling(density);
	arma::mat B = (1.0 + BExternalGlobaltemp) - 2.0 * abs2ToFilling(density);

	
	arma::mat A_x=0.0*abs2ToFilling(density);
	A_x.zeros();
	arma::mat A_y=A_x;	

	auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
	auto y_origin = round(my);
	for (iint i = 0; i < Nx; i++) {
		A_x(i, y_origin) = 0.0;
	}
	
    for (iint j = y_origin + 1; j < Ny; j++) {
		for (iint i = 0; i < Nx; i++) {
			A_x(i, j) = A_x(i, j-1) + B(i, j-1) * ay;
		}
	}
	
    for (iint j = y_origin - 1; j >= 0; j--) {
		for (iint i = 0; i < Nx; i++) {
			A_x(i, j) = A_x(i, j+1) - B(i, j) * ay;
		}
	}
	
	//cout << "debug: " << accu(A_x) << endl;

	A_x_mat = A_x;
	
	A_y_mat = A_y;
	
}


//to do: use 4 points average
void calcA_Landau_x_refined(const arma::mat &density /*wave function squared*/, const double& BExternalGlobaltemp) {
	using namespace VECTORPOT;//yayun, 240101, sure necessary to reuse but not sure how to avoid reuse
	
	//arma::mat B_reduced = (1.0 + BExternalGlobaltemp) -2.0*abs2ToFilling(density);
	arma::mat B = (1.0 + BExternalGlobaltemp) - 2.0 * abs2ToFilling(density);

	
	arma::mat A_x = 0.0 * abs2ToFilling(density);
	A_x.zeros();
	arma::mat A_y=A_x;	

	auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
	auto y_origin = round(my);
	for (iint i = 0; i < Nx; i++) {
		A_x(i, y_origin) = 0.0;
	}
	
    for (iint j = y_origin + 1; j < Ny; j++) {
		for (iint i = 0; i < Nx; i++) {
			A_x(i, j) = A_x(i, j-1) - (B(i, j-1) + B(i, j)) / 2.0 * ay;
		}
	}
	
    for (iint j = y_origin - 1; j >= 0; j--) {
		for (iint i = 0; i < Nx; i++) {
			A_x(i, j) = A_x(i, j+1) + (B(i, j) + B(i, j+1)) / 2.0 * ay;
		}
	}
	
	//cout << "accu(A_x) debug: " << accu(A_x) << endl;
	//cout << "accu(A_y) debug: " << accu(A_y) << endl;

	A_x_mat = A_x;
	
	A_y_mat = A_y;
	
}


/* for debugging purpose, add the contribution from the circumcircle of the square lattice patch
 * to be run after each calcA() call
 */
void reduntA() {
	using namespace arma;
	using namespace VECTORPOT;

	mat gx = A_x_mat;
	mat gy = A_y_mat;

	mat zd(Nx,Ny);
	zd.fill(0.0);
	calcA(zd, BExternalGlobal);
	mat zx = A_x_mat;
	mat zy = A_y_mat;

	A_x_mat = gx - zx;
	A_y_mat = gy - zy;

	mat bx(Nx,Ny);
	mat by(Nx,Ny);

	bx.fill(0.0);
	by.fill(0.0);
	
	const auto c = 1.0/2;

	auto mx=(Nx-1)/2, my=(Ny-1)/2;//position index of the middle point in the x and y directions

	for(iint j=0; j<Ny; j++) {
		for(iint i=0; i<Nx; i++) {
			double r = sqrt((i-mx)*(i-mx)*ax*ax + (j-my)*(j-my)*ay*ay);
			
			bx(i,j)=c*(my-j)*ay;//analytic solution for uniform magnetic field in regions of the rectangle plus the excluded circumference circle
			by(i,j)=c*(i-mx)*ax;
		}
	}
	A_x_mat += bx;
	A_y_mat += by;

	bx.save("bx-1.csv", arma::csv_ascii);
	by.save("by-1.csv", arma::csv_ascii);

	zx.save("zx-1.csv", arma::csv_ascii);
	zy.save("zy-1.csv", arma::csv_ascii);
	gx.save("gx-1.csv", arma::csv_ascii);
	gy.save("gy-1.csv", arma::csv_ascii);
	
	gx = gx-bx;
	gy = gy-by;
}


/* create an extra vector potential from the contribution of the circumcircle of the rectangle lattice patch; save it for once for use externally by matlab, 
but do not need to load it because it can be calculated fresh (if modified during a parameter) for each new table parameter of the external magnetic field; */
void AddCircumcircleA(const bool& circleA, const double & BExternalGlobaltemp, const int whichexternaltemp, const double ionFillingtemp) {//be careful that A_x_mat will change whenever calcA() is called;
	if (circleA==true){
	    using namespace arma;
	    using namespace VECTORPOT;
		
		if (whichexternaltemp==0){
			//save A from only effective magnetic field on the rectangle to gx and gy (i.e. A from the real self-consistent system which is now under the iterative calculation)
			mat gx = A_x_mat;
			mat gy = A_y_mat;
			//calculate A from only real uniform magnetic field on the rectangle and save to zx and zy
			mat zd(Nx, Ny);
			zd.fill(0.0);//meaning filling factor 0.0, such that B=1 is external
			calcA(zd, 0.0);//do not consider BExternalGlobaltemp here, taken to be zero; last step will consider extra magnetic field
			mat zx = A_x_mat;
			mat zy = A_y_mat;
			//analytic solution of A from real uniform magnetic field on the circle (including the rectangle) and save to bx and by
			mat bx(Nx, Ny);
			mat by(Nx, Ny);
			bx.fill(0.0);
			by.fill(0.0);

			const auto c = 1.0 / 2;
			auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
			for (iint j = 0; j < Ny; j++) {
				for (iint i = 0; i < Nx; i++) {
					double r = sqrt((i - mx) * (i - mx) * ax * ax + (j - my) * (j - my) * ay * ay);

					bx(i, j) = c * (my - j) * ay;//analytic solution for uniform magnetic field in regions of the rectangle plus the excluded circumference circle
					by(i, j) = c * (i - mx) * ax;
				}
			}
			double coeff = 0.0;
			coeff = 1.0 * (1.0 + BExternalGlobaltemp);			
			//the external circumcircle A: bx - zx (i.e. analytic minus rectangular contribution)
			bx = (bx - zx) * coeff;
			by = (by - zy) * coeff;
			//sum up rectange A: gx in the selfconsistent system and the external circumcircle A
			A_x_mat = gx + bx;
			A_y_mat = gy + by;			
		}else if(whichexternaltemp==1){
			//save A from only effective magnetic field on the rectangle to gx and gy (i.e. A from the real self-consistent system which is now under the iterative calculation)
			mat gx = A_x_mat;
			mat gy = A_y_mat;
			//calculate A from only real uniform magnetic field on the rectangle and save to zx and zy
			mat zd(Nx, Ny);
			zd.fill(0.0);//meaning filling factor 0.0, such that B=1 is external
			calcA(zd, 0.0);//do not consider BExternalGlobaltemp here, taken to be zero; last step will consider extra magnetic field
			mat zx = A_x_mat;
			mat zy = A_y_mat;
			//analytic solution of A from real uniform magnetic field on the circle (including the rectangle) and save to bx and by
			mat bx(Nx, Ny);
			mat by(Nx, Ny);
			bx.fill(0.0);
			by.fill(0.0);

			const auto c = 1.0 / 2;
			auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
			for (iint j = 0; j < Ny; j++) {
				for (iint i = 0; i < Nx; i++) {
					double r = sqrt((i - mx) * (i - mx) * ax * ax + (j - my) * (j - my) * ay * ay);

					bx(i, j) = c * (my - j) * ay;//analytic solution for uniform magnetic field in regions of the rectangle plus the excluded circumference circle
					by(i, j) = c * (i - mx) * ax;
				}
			}
			double coeff = 0.0;
			coeff = (1.0 - 2.0 * ionFillingtemp) * (1.0 + BExternalGlobaltemp);
			//the external circumcircle A: bx - zx (i.e. analytic minus rectangular contribution)
			bx = (bx - zx) * coeff;
			by = (by - zy) * coeff;
			//sum up rectange A: gx in the selfconsistent system and the external circumcircle A
			A_x_mat = gx + bx;
			A_y_mat = gy + by;
		}else if(whichexternaltemp == 2){
			//save A from only effective magnetic field on the rectangle to gx and gy (i.e. A from the real self-consistent system which is now under the iterative calculation)
			mat gx = A_x_mat;
			mat gy = A_y_mat;
			//calculate A from only real uniform magnetic field on the rectangle and save to zx and zy
			mat zd(Nx, Ny);
			zd.fill(0.0);//meaning filling factor 0.0, such that B=1 is external
			calcA(zd, 0.0);//do not consider BExternalGlobaltemp here, taken to be zero; last step will consider extra magnetic field
			mat zx = A_x_mat;
			mat zy = A_y_mat;
			//analytic solution of A from real uniform magnetic field on the infinite long rectangle along the x direction of width Ly (including the rectangle) and save to bx and by
			mat bx(Nx, Ny);
			mat by(Nx, Ny);
			bx.fill(0.0);
			by.fill(0.0);

			const auto c = 1.0;// this is now not equal to 1/2
			auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
			for (iint j = 0; j < Ny; j++) {
				for (iint i = 0; i < Nx; i++) {
					double r = sqrt((i - mx) * (i - mx) * ax * ax + (j - my) * (j - my) * ay * ay);

					bx(i, j) = c * (my - j) * ay;//analytic solution for uniform magnetic field in regions of the rectangle plus the excluded two semi-infinite leads
					by(i, j) = 0;//c * (i - mx) * ax;
				}
			}
			double coeff = 0.0;
			coeff = (1.0 - 2.0 * ionFillingtemp) * (1.0 + BExternalGlobaltemp);
			//the external circumcircle A: bx - zx (i.e. analytic minus rectangular contribution)
			bx = (bx - zx) * coeff;
			by = (by - zy) * coeff;
			//sum up rectange A: gx in the selfconsistent system and the external circumcircle A
			A_x_mat = gx + bx;
			A_y_mat = gy + by;
		}else if(whichexternaltemp == 3){
			//save A from only effective magnetic field on the rectangle to gx and gy (i.e. A from the real self-consistent system which is now under the iterative calculation)
			mat gx = A_x_mat;
			mat gy = A_y_mat;
			//calculate A from only real uniform magnetic field on the rectangle and save to zx and zy
			mat zd(Nx, Ny);
			zd.fill(0.0);//meaning filling factor 0.0, such that B=1 is external
			calcA(zd, 0.0);//do not consider BExternalGlobaltemp here, taken to be zero; last step will consider extra magnetic field
			mat zx = A_x_mat;
			mat zy = A_y_mat;
			//analytic solution of A from real uniform magnetic field on the infinite long rectangle along the x direction of width Ly (including the rectangle) and save to bx and by
			mat bx(Nx, Ny);
			mat by(Nx, Ny);
			bx.fill(0.0);
			by.fill(0.0);

			const auto c = 1.0;// this is now not equal to 1/2
			auto mx = (Nx - 1) / 2, my = (Ny - 1) / 2;//position index of the middle point in the x and y directions
			for (iint j = 0; j < Ny; j++) {
				for (iint i = 0; i < Nx; i++) {
					double r = sqrt((i - mx) * (i - mx) * ax * ax + (j - my) * (j - my) * ay * ay);

					bx(i, j) = 0;//c * (my - j) * ay;//analytic solution for uniform magnetic field in regions of the rectangle plus the excluded two semi-infinite leads
					by(i, j) = c * (i - mx) * ax;
				}
			}
			double coeff = 0.0;
			coeff = (1.0 - 2.0 * ionFillingtemp) * (1.0 + BExternalGlobaltemp);
			//the external circumcircle A: bx - zx (i.e. analytic minus rectangular contribution)
			bx = (bx - zx) * coeff;
			by = (by - zy) * coeff;
			//sum up rectange A: gx in the selfconsistent system and the external circumcircle A
			A_x_mat = gx + bx;
			A_y_mat = gy + by;
		}


	}
}

/*****************************************
 * Hartree potential, i.e. Coulomb potenial.
 *****************************************/
namespace COULOMB {
	static arma::mat greenCoulomb; // Green's functions for 1/r.
	static arma::mat VCoulomb; // Coulomb potential, the convolution result.	
	/*
	 * Electrostatic energy of a charged 1x1 square is
	 * 2/3 (2-2 Sqrt[2]+3 ArcSinh[1]+3 Log[1+Sqrt[2]])
	 * (* Great. Here are some rectangles... *)
	 */
	const double coulombSquare = 2.0/3*(2-2*sqrt(2)+3*asinh(1)+3*log(1+sqrt(2)))/sqrt(ax*ay);
	double softCoulomb = coulombSquare;
	
	/* MKL_VSL */
	VSLConvTaskPtr taskCoulomb;


};

/* Fill in Green's functions */
void initGreenCoulomb() {
	using namespace COULOMB;
	using arma::mat;

	/* filling Green function matrices */
	int ayl = 2*Ny-1;
	int axl = 2*Nx-1;
	int mx=(axl-1)/2, my=(ayl-1)/2;
	greenCoulomb = mat(axl,ayl);
	VCoulomb = mat(Nx,Ny);
	
	for(iint i=0; i<axl; i++) {
		for(iint j=0; j<ayl; j++) {
			double r2 = (i-mx)*(i-mx)*ax*ax*DielecX_vs_DielecY + (j-my)*(j-my)*ay*ay;
			greenCoulomb(i,j)=sqrt(1/r2);
		}
	}
	greenCoulomb(mx,my) = softCoulomb;// The potential at the point center, created by the charge at the point center
}

/* MKL convolution initialize for Coulomb */
void createConvTaskCoulomb() {
	using namespace COULOMB;
	auto mode = VSL_CONV_MODE_AUTO;
	auto green = reinterpret_cast<const double*>(greenCoulomb.memptr());
	MKL_INT dims = 2;
	MKL_INT xshape[2] = {(MKL_INT)greenCoulomb.n_rows,(MKL_INT)greenCoulomb.n_cols};
	MKL_INT yshape[2] = {(MKL_INT)VCoulomb.n_rows,(MKL_INT)VCoulomb.n_cols};
	MKL_INT xstride[2] = {1,xshape[0]};
	MKL_INT start[2] = {(MKL_INT)(greenCoulomb.n_rows/2),(MKL_INT)(greenCoulomb.n_cols/2)};

	auto status = vsldConvNewTaskX(&taskCoulomb, mode, dims, xshape, yshape, yshape, green, xstride);
	vslerror(status, "Failed vsldConvNewTaskX for taskCoulomb");
	status = vslConvSetStart(taskCoulomb, start);
	vslerror(status, "Failed vsldConvSetStart for taskCoulomb");
}

/* Convolve to get Coulomb potential */
void calcCoulomb(const arma::mat &density) {
	using namespace COULOMB;

	auto vc = reinterpret_cast<double*>(VCoulomb.memptr());
	auto rho = reinterpret_cast<const double*>(density.memptr());
	MKL_INT ystride[2] = {1,(MKL_INT)density.n_rows};

	auto status = vsldConvExecX(taskCoulomb,rho,ystride,vc,ystride);
	vslerror(status, "Failed vsldConvExecX for VCoulomb");
}








/*****************************************
 * the V_T term as two parts
 *****************************************/
namespace VTSTAR {
	arma::mat greenVT_x, greenVT_y; // Green's functions.
	arma::mat VT_x_mat, VT_y_mat; // VTx & VTy. Convolution results.
	/* MKL_VSL */
	VSLConvTaskPtr taskVT_x, taskVT_y;
	const auto VTc = 0.001*-4*Dielec/alpha_m;
};


void initGreenVT() {
	using namespace VTSTAR;

	/* filling Green function matrices */
	auto ayl = 2*Ny-1;
	auto axl = 2*Nx-1;
	auto mx=(axl-1)/2, my=(ayl-1)/2;
	greenVT_x = arma::mat(axl,ayl);
	greenVT_y = arma::mat(axl,ayl);
	
	const auto c = ax*ay;// Why there is a c=ax*ay here?

	for(iint j=0; j<ayl; j++) {
		for(iint i=0; i<axl; i++) {
			double r2 = (i-mx)*(i-mx)*ax*ax + (j-my)*(j-my)*ay*ay;
			
			greenVT_x(i,j)=-VTc*c*(my-j)*ay/r2;
			greenVT_y(i,j)=-VTc*c*(i-mx)*ax/r2;
		}
	}
	greenVT_x(mx,my) = 0;
	greenVT_y(mx,my) = 0;
	
	VT_x_mat = arma::mat(Nx,Ny);
	VT_y_mat = arma::mat(Nx,Ny);
}

void createConvTaskVT() {
	using namespace VTSTAR;

	auto	mode	   = VSL_CONV_MODE_AUTO;
	auto	gVT_x	   = reinterpret_cast<const double*>(greenVT_x.memptr());
	auto	gVT_y	   = reinterpret_cast<const double*>(greenVT_y.memptr());
	MKL_INT dims	   = 2;
	MKL_INT xshape[2]  = {(MKL_INT)greenVT_x.n_rows,(MKL_INT)greenVT_x.n_cols};
	MKL_INT yshape[2]  = {(MKL_INT)VT_x_mat.n_rows,(MKL_INT)VT_x_mat.n_cols};
	MKL_INT xstride[2] = {1,xshape[0]};
	/* The default start is {0,0,...}, 0-based. */
	MKL_INT start[2]   = {(MKL_INT)(greenVT_x.n_rows/2),(MKL_INT)(greenVT_x.n_cols/2)};

	auto status = vsldConvNewTaskX(&taskVT_x, mode, dims, xshape, yshape, yshape, gVT_x, xstride);
	vslerror(status, "Failed vsldConvNewTaskX for VT_x");
	status = vslConvSetStart(taskVT_x, start);
	vslerror(status, "Failed vslConvSetStart for VT_x");
	
	status = vsldConvNewTaskX(&taskVT_y, mode, dims, xshape, yshape, yshape, gVT_y, xstride);
	vslerror(status, "Failed vsldConvNewTaskX for VT_y");
	status = vslConvSetStart(taskVT_y, start);
	vslerror(status, "Failed vslConvSetStart for VT_y");
}


/* Convolve to get VT */
void calcVT(const arma::mat &px_mat, const arma::mat &py_mat /* p-qA, velocity matrix */) {
	using namespace VTSTAR;

	MKL_INT ystride[2] = {1,static_cast<int>(px_mat.n_rows)};
	
	auto px = reinterpret_cast<const double*>(px_mat.memptr());
	auto py = reinterpret_cast<const double*>(py_mat.memptr());
	auto VT_x = reinterpret_cast<double*>(VT_x_mat.memptr());
	auto VT_y = reinterpret_cast<double*>(VT_y_mat.memptr());

	auto status = vsldConvExecX(taskVT_x,px,ystride,VT_x,ystride);
	vslerror(status, "Failed vsldConvExecX for VT_x");

	status = vsldConvExecX(taskVT_y,py,ystride,VT_y,ystride);
	vslerror(status, "Failed vsldConvExecX for VT_y");
}













#endif
