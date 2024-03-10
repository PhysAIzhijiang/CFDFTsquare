// This file is under public domain.
#ifndef SYSPARA_HPP
#define SYSPARA_HPP

#include "commoninclude.hpp"

enum VERBOSITY{ Quiet, Normal, Info, Debug, DumpAll };
const enum VERBOSITY verbosity = Debug;
const std::string prefix = ""; // prefix to file names
//const int pfrequency = 500, eesfrequency = 500, evsfrequency = 2000, denfrequency = 500;
const int pfrequency = 2000, eesfrequency = 50, evsfrequency = 2000, denfrequency = 50;
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

enum SHAPE {Circle, Rectangle, FullRectangle};
const enum SHAPE shape = FullRectangle;//Circle; 
const int shape_fill_num=2;
//changes 20240116: ionfilling, calcA_Landau_x: 0.0*abs2ToFilling(density) 
const double ionFilling = 2.0/5.0;//1.0/3; // initial filling factor, over the positive ion region.
//arma::rowvec BExternalGlobalSet = arma::linspace<arma::rowvec>(0.0,-0.06,400);
arma::rowvec BExternalGlobalSet = arma::linspace<arma::rowvec>(0.0,0.05,400);
//arma::rowvec BExternalGlobalSet = arma::linspace<arma::rowvec>(0.0,0.03,200);
//arma::linspace<arma::rowvec>(0,0.04,161);
//{ 0.0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.01};//arma::linspace<rowvec>(0, 0.02, 41);
//{ 0.0, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, 0.0475, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0625, 0.065, 0.0675, 0.07, 0.0725, 0.075, 0.0775, 0.08, 0.0825, 0.085, 0.0875, 0.09, 0.0925, 0.095, 0.0975, 0.10};
//{ 0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.20, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, 0.235, 0.24 };
//{ 0.0, -0.005, -0.01, -0.015, -0.02, -0.025, -0.03, -0.035, -0.04, -0.045, -0.05, -0.055, -0.06, -0.065, -0.07, -0.075, -0.08, -0.085, -0.09, -0.095, -0.10, -0.105, -0.11, -0.115, -0.12, -0.125, -0.13, -0.135, -0.14, -0.145, -0.15, -0.155, -0.16, -0.165, -0.17, -0.175, -0.18, -0.185, -0.19, -0.195, -0.20, -0.205, -0.21, -0.215, -0.22, -0.225, -0.23, -0.235, -0.24};
//
//BExternalGlobalSet = BExternalGlobalSet * (-1.0);
//{ 0.0, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24 };
//{ 0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24 };

//{ 0.0, -0.03, -0.06, -0.09, -0.12, -0.15, -0.18, -0.21, -0.24 };//{ 0.0, -0.02, -0.04, -0.06, -0.08, -0.1, -0.12, -0.14, -0.16 };//{ 0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16 };//, 0.08, 0.09, 0.10, 0.12, 0.14, 0.16
double BExternalGlobal = 0.0;
arma::rowvec MyTemperatureSet = {0.000001};//{0.001};//{0.001, 0.0001, 0.00001};//{0.005, 0.001, 0.0001, 0.00001}; //{ 0.05, 0.01, 0.005, 0.001 };
double kBT = 0.001;
arma::rowvec MySecondSet = BExternalGlobalSet;
int TemperatureNum = MyTemperatureSet.n_elem; 
int MySecondSetNum = MySecondSet.n_elem;
int KeepUsingInitialConditions = 1;
int LoadAndReiterateConvergedDada = 0;
/*
 * Knobs for different factors of the functional
 */
const double scaleHartree = 0.1;//0.2;
const double disorderStrength = 0.03;
const double scaleXC = 0.0; //1.0
const double scaleIon = scaleHartree;
const double scaleVT = 0.0;
const double gradcoe= 0.0;
//const int interceptXC = 2.0;
/*****************************************
 * computational parameters
 *****************************************/

/* lattice parameters *//* system size, in magnetic length l_B*/// resolution per l_B
//const int Lx = 6, Ly = 3, Nrx = 2, Nry = 2;
//const int Lx = 10, Ly = 10, Nrx = 4, Nry = 4;
//const int Lx = 30, Ly = 15, Nrx = 2, Nry = 2;
//const int Lx = 15, Ly = 15, Nrx = 4, Nry = 4;
//const int Lx = 10, Ly = 100, Nrx = 10, Nry = 10;
//const int Lx = 40, Ly = 30, Nrx = 10, Nry = 10;
//const int Lx = 40, Ly = 30, Nrx = 3, Nry = 3;
//const int Lx = 60, Ly = 45, Nrx = 2, Nry = 2;
//const int Lx = 90, Ly = 30, Nrx = 2, Nry = 2;
//const int Lx = 120, Ly = 30, Nrx = 2, Nry = 2;
//const int Lx = 140, Ly = 40, Nrx = 2, Nry = 2;
//const int Lx = 125, Ly = 40, Nrx = 2, Nry = 2;
//const int Lx = 200, Ly = 64, Nrx = 2, Nry = 2;
//const int Lx = 150, Ly = 60, Nrx = 2, Nry = 2;
//const int Lx = 120, Ly = 60, Nrx = 2, Nry = 2;
//const int Lx = 100, Ly = 70, Nrx = 2, Nry = 2;
//const int Lx = 100, Ly = 70, Nrx = 1, Nry = 1;
const int Lx = 200, Ly = 140, Nrx = 1, Nry = 1;
//const int Lx = 120, Ly = 80, Nrx = 2, Nry = 2;
//const int Lx = 120, Ly = 100, Nrx = 2, Nry = 2;
//const int Lx = 100, Ly = 50, Nrx = 2, Nry = 2;
//const int Lx = 40, Ly = 30, Nrx = 15, Nry = 15;
//const int Lx = 30, Ly = 30, Nrx = 6, Nry = 6;
//const int Lx = 35, Ly = 35, Nrx = 4, Nry = 4;
//const int Lx = 35, Ly = 35, Nrx = 6, Nry = 6;
//const int Lx = 35, Ly = 35, Nrx = 4, Nry = 4;
//const int Lx = 35, Ly = 35, Nrx = 10, Nry = 10;
//const int Lx = 120, Ly = 45, Nrx = 3, Nry = 3;
//const int Lx = 120, Ly = 45, Nrx = 3, Nry = 3;
//const int Lx = 88, Ly = 33, Nrx = 3, Nry = 3;
//const int Lx = 45, Ly = 45, Nrx = 4, Nry = 4;
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
const auto CyclotronGap = 2 * KE_scale * Bu_vs_B / (sqrt(Mx_vs_My)); // omega_c, CyclotronGap is heavily used in guessing eigenenergy bounds. But Bu is not used.
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
bool knowFeastPara = false;//the quick way to set initial guess of energy range.
double knowEnergylow = -0.132909;
double knowEnergyhigh = 0.216799;
int knowNsubspace = 2500;
bool useGuessSubspaceZeroth = false;
/*
Temperature related settings. The parameters are independent of the values of kBT and so on, so do not need updates
*/
bool isThereTemperature = true;
const double deltaElectronNumToleranceReal = 0.0000000001; //Try to find the most accurate occu possible, but also accept if a flexible bar is met
const double deltaElectronNumTolerance = 0.0000000000000000000001;// This is to reach a meaninglessly small numerical error
const double percentDiffRatio = 0.0;
//double kBT = 0.001;
double occuMinimumFirst = 0.001;//It seems that 10^(-15) is numerical instability, so use it as my Cmin in imagination
double occuMinimumEkBTFirst = log(1.0/ occuMinimumFirst -1.0);//=6.9068; EF + occuMinimumEkBTFirst * kBt gives this occupation;
double occuMinimumSecond = occuMinimumFirst *0.01;
double occuMinimumEkBTSecond = log(1.0 / occuMinimumSecond - 1.0);//=11.5129; EF + occuMinimumEkBTSecond * kBt gives this occupation;
/*
 * Iteration parameters
 */
const double Temperature = 0.0; // initial temperature for annealing
double UpdateRatio = 0.05; //wavefunction density update ratio; this number here is not used, because it is erased/assigned in each for loop of table parameter
double UpdateRatioVT = UpdateRatio; // amount of VT potential to be updated
const int Niter = 1000, lookBackStepNum = 10;
const double UpdateRatioBase = UpdateRatio; // amount of total wavefunction density to be updated
const double UpdateRatioLowBound = 0.01;
const double UpdateRatioUpperBound = 0.2;
//const double UpdateRatioSmall = 0.01;
double UpdateRatioTemp = 0.01;
/*
 * Stopping criteria
 */
const double enterBigThreshold = 0.1;
const double relDensityThreshold = 0.0001;
const double BigNumThreshold = 0.1;
const double chemicalNumThreshold = 0.01;
const double KineticThreshold = 0.01; 
int OccuEigenStorage = 100;//column size of TableOccuEigen
int HistoryDetailStorage = 3000;//column size of TableHistoryDetail
/*
*effects of external environment
*/
bool externalPotentialGivenFromOutside = false;
bool useCircumcircleA = false; //whether use an extra vector potential from the contribution of the circumcircle of the rectangle lattice patch;
const int whichexternal = 2;//external useCircumcircleA contribution: 0 for real B and, 1 for uniform B* from perfect filling of ionfilling; 2 for infinite long 
bool EfModeFlag = false;//true;  Yayun 240115, maybe false is more adiabatic?
bool usePer = false;//whether use periodic boundary conditions along x direction

#endif
