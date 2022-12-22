// This file is under public domain.
#ifndef SYSPARA_HPP
#define SYSPARA_HPP

#include "commoninclude.hpp"

enum VERBOSITY{ Quiet, Normal, Info, Debug, DumpAll };
const enum VERBOSITY verbosity = Debug;
const std::string prefix = ""; // prefix to file names
//const int pfrequency = 500, eesfrequency = 500, evsfrequency = 2000, denfrequency = 500;
const int pfrequency = 2000, eesfrequency = 2000, evsfrequency = 2000, denfrequency = 200;
//const int pfrequency = 100, eesfrequency = 50, evsfrequency = 100, denfrequency = 20;
//const int pfrequency = 30, eesfrequency = 30, evsfrequency = 1000, denfrequency = 30;
//const int pfrequency = 1000, eesfrequency = 1000, evsfrequency = 1000, denfrequency = 2;
//const int pfrequency = 2, eesfrequency = 2, evsfrequency = 2, denfrequency = 2;
/*save cycles*/
const int save_cycleNum1 = 2/*density*/, save_cycleNum2 = 2/*potentials*/, save_cycleNum3 = 3/*not used*/;
/*****************************************
 * Physical parameters
 * The energy scale is e^2/l_B/epsilon.
 * The length scale is l_B.
 *****************************************/

const double alpha_m = 0.08*1.0;	// effective mass
const double Mx_vs_My = 1.0; // Anisotropy; "vs" stands for "over".
const double Dielec = 13.6; // GaAs
const double DielecX_vs_DielecY = 1.0;
const auto KE_scaleReferenceB = 0.001031 * Dielec / alpha_m; //coefficient to convert kinteic energy to e^2/l_B/epsilon
double KE_scale = 0.001031 * Dielec / alpha_m; //coefficient to convert kinteic energy to e^2/l_B/epsilon

enum SHAPE {Circle, Rectangle};
const enum SHAPE shape = Rectangle;//Circle; 

const double ionFilling = 1.0/3; // initial filling factor, over the positive ion region.
arma::rowvec BExternalGlobalSet = { 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24 };
//{ 0.0, -0.03, -0.06, -0.09, -0.12, -0.15, -0.18, -0.21, -0.24 };//{ 0.0, -0.02, -0.04, -0.06, -0.08, -0.1, -0.12, -0.14, -0.16 };//{ 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16 };//, 0.08, 0.09, 0.10, 0.12, 0.14, 0.16
double BExternalGlobal = 0.0;
arma::rowvec MyTemperatureSet = { 0.01, 0.005, 0.001, 0.0001, 0.00001}; //{ 0.05, 0.01, 0.005, 0.001 };
double kBT = 0.001;
arma::rowvec MySecondSet = BExternalGlobalSet;
int TemperatureNum = MyTemperatureSet.n_elem; 
int MySecondSetNum = MySecondSet.n_elem;
int KeepUsingInitialConditions = 0;
int LoadAndReiterateConvergedDada = 0;
/*
 * Knobs for different factors of the functional
 */
const double scaleHartree = 1.0;//0.2;
const double disorderStrength = 0.00;
const double scaleXC = 1.0; //1.0
const double scaleIon = scaleHartree;
const double scaleVT = 1.0;
const double gradcoe= 0.0;
//const int interceptXC = 2.0;
/*****************************************
 * computational parameters
 *****************************************/

 /* lattice parameters *//* system size, in magnetic length l_B*/// resolution per l_B
//const int Lx = 10, Ly = 10, Nrx = 3, Nry = 3;
//const int Lx = 15, Ly = 15, Nrx = 3, Nry = 3;
//const int Lx = 15, Ly = 15, Nrx = 4, Nry = 4;
//const int Lx = 30, Ly = 30, Nrx = 3, Nry = 3;
//const int Lx = 30, Ly = 30, Nrx = 6, Nry = 6;
//const int Lx = 35, Ly = 35, Nrx = 4, Nry = 4;
//const int Lx = 35, Ly = 35, Nrx = 6, Nry = 6;
//const int Lx = 35, Ly = 35, Nrx = 4, Nry = 4;
//const int Lx = 35, Ly = 35, Nrx = 10, Nry = 10;
//const int Lx = 35, Ly = 35, Nrx = 4, Nry = 4;
const int Lx = 45, Ly = 45, Nrx = 4, Nry = 4;
//const int Lx = 70, Ly = 70, Nrx = 4, Nry = 4;
//const int Lx = 70, Ly = 70, Nrx = 3, Nry = 3;
//const int Lx = 70, Ly = 70, Nrx = 5, Nry = 5;
//const int Lx = 70, Ly = 70, Nrx = 6, Nry = 6;
//const int Lx = 35, Ly = 35, Nrx = 6, Nry = 6;
//const int Lx = 30, Ly = 30, Nrx = 3, Nry = 3;
//const int Lx = 30, Ly = 30, Nrx = 6, Nry = 6;
//const int Lx = 30, Ly = 30, Nrx = 10, Nry = 10;
//const int Lx = 30, Ly = 30, Nrx = 15, Nry = 15;
//const int Lx = 50, Ly = 50, Nrx = 6, Nry = 6; 
const iint Nx = Lx * Nrx, Ny = Ly * Nry, Nxy = Nx * Ny; // number of array points
iint Np = 60;//55; // May be overriden if continuing from old SeaFilling files.
       

/* grid cell length in units of magnetic length */
const auto ax = 1.0 / Nrx, ay = 1.0 / Nry;

/*Leave the issues of Bu_vs_B here for now*/
const auto Bu_vs_B = (1 - 2 * ionFilling)+ 0.0;//BExternalGlobal // Bu/B , vs for division here. Yayun,0502, this Bu is not used, this quantity is not necessary
const auto CyclotronGap = 2 * KE_scale * Bu_vs_B / (sqrt(Mx_vs_My)); // omega_c, CyclotronGap is heavily used in guessing eigenenergy bounds
const auto Nrxu = Nrx / sqrt(std::abs(Bu_vs_B)), Nryu = Nry / sqrt(std::abs(Bu_vs_B)); // resolution per l_Bu; There is no B_u used in this code.
const auto Nstate_per_lambda = Bu_vs_B*Lx*Ly/2/Pi;

/*
 * common helper functons
 */
const auto filling_vs_abs2 = (2*Pi*Nrx*Nry);

template<typename T>
T abs2ToFilling(const T wfsq /*sum of wave function squared*/) {
	T a = filling_vs_abs2 * wfsq;
	return a;
}
template<typename T>
T fillingToAbs2(const T filling) {
	T a = filling / filling_vs_abs2;
	return a;
}
template<typename T>
T abs2ToFillingRealB(const T wfsq /*sum of wave function squared*/, const double deltaB) {
	const auto filling_vs_abs2RealB = (2 * Pi * Nrx * Nry)/(1.0 + deltaB);
	T a = filling_vs_abs2RealB * wfsq;
	return a;
}
template<typename T>
T fillingToAbs2RealB(const T filling, const double deltaB) {
	const auto filling_vs_abs2RealB = (2 * Pi * Nrx * Nry) / (1.0 + deltaB);
	T a = filling / filling_vs_abs2RealB;
	return a;
}

template<typename T>
T multiply_coefficient(const T matrix, double coefficient) {
	T a = coefficient * matrix;
	return a;
}//not used yet

/*
 * Additional extra potential.
 */
#include "extraV.hpp"
const double confinementHeight = 0.0;
const auto addConfinement = [](MMSpmatR_c64 &h, double v){diskWell(h,v);}; // function to generate confinement potential
const std::string confinementName = "diskWell";
bool useAnalyticBackgroundPotential = true;
const double insideDeltaBarrierNp = Np/4.0;
/*
 * Initial state parameters.
 * In a continued run, these parameters are not used.
 *
 * Electron number is set by ion radius, since system is made neutral.
 * Later electron number is made integer, at the sacrifice of filling accuracy.
 */
const double geometryYoverX = 1.0;
auto ionRadius = sqrt(Np / fillingToAbs2(ionFilling) / Pi); // Radius of disk region for the initial state and the ion background.
auto xLen = sqrt(Np / fillingToAbs2(ionFilling) / geometryYoverX); // length of side x region for the initial state and the ion background.
auto yLen = xLen * geometryYoverX;
const double initElectronFilling = ionFilling;
/** end initial state parameters **/

/*
 * State or solution parameters
 */

/* Hard bounds for eigensubspace size during initial diagonalization. */
const iint NsubspaceHardMax = 1<<20;
const iint NsubspaceHardMin = 1<<8;
/* initial bounds for eigenvalues */
double Energylow = -10.0;
double Energyhigh = 1.0;
bool useBoundFromPotential = true;
bool knowFeastPara = false;
double knowEnergylow = -0.2;
double knowEnergyhigh = 0.1;
int knowNsubspace = 50;
/*
Temperature related settings. The parameters are independent of the values of kBT and so on, so do not need updates
*/
bool isThereTemperature = true;
const double deltaElectronNumToleranceReal = 0.0000000001; //Try to find the most accurate occu possible, but also accept if a flexible bar is met
const double deltaElectronNumTolerance = 0.0000000000000000000001;// This is to reach a meaninglessly small numerical error
const double percentDiffRatio = 0.0;
//double kBT = 0.001;
double occuMinimumFirst = 0.001;//It seems that 10^(-15) is numerical instability, so use it as my Cmin in imagination
double occuMinimumEkBTFirst = log(1.0/ occuMinimumFirst -1.0);//EF + occuMinimumEkBTFirst * kBt gives this occupation;
double occuMinimumSecond = occuMinimumFirst *0.01;
double occuMinimumEkBTSecond = log(1.0 / occuMinimumSecond - 1.0);// EF + occuMinimumEkBTSecond * kBt gives this occupation;
/*
 * Iteration parameters
 */
const double Temperature = 0.0; // initial temperature for annealing
double UpdateRatio = 0.05; //wavefunction density update ratio; this number here is not used, because it is erased/assigned in each for loop of table parameter
double UpdateRatioVT = UpdateRatio; // amount of VT potential to be updated
const int Niter = 2001, lookBackStepNum = 10;
const double UpdateRatioBase = UpdateRatio; // amount of total wavefunction density to be updated
const double UpdateRatioLowBound = 0.01;
const double UpdateRatioUpperBound = 0.2;
/*
 * Stopping criteria
 */
const double enterBigThreshold = 0.1;
const double relDensityThreshold = 0.000001;
const double BigNumThreshold = 0.1;
const double chemicalNumThreshold = 0.01;
const double KineticThreshold = 0.01; 
int OccuEigenStorage = 100;//column size of TableOccuEigen
int HistoryDetailStorage = 3000;//column size of TableHistoryDetail

bool externalPotentialGivenFromOutside = false;
bool useCircumcircleA = true; 
bool EfModeFlag = true;

#endif
