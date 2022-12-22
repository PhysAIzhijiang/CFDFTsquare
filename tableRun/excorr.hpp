// This file is under public domain.
#ifndef EXCORR_HPP
#define EXCORR_HPP

#include "commoninclude.hpp"
#include "syspara.hpp"
//#include "lattice.h"

/*
 * Exchange correlation.
 */
double VxcRef(const double);
double ExcRef(const double);

double VxcPaper(const double);
double ExcPaper(const double);


double (*Vxc)(const double) = VxcPaper;
double (*Exc)(const double) = ExcPaper;
//double (*Vxc)(const double) = VxcRef;
//double (*Exc)(const double) = ExcRef;

double VxcPaperRealFilling(const double);
double ExcPaperRealFilling(const double, const double);

double (*VxcRealFilling)(const double) = VxcPaperRealFilling;
double (*ExcRealFilling)(const double, const double) = ExcPaperRealFilling;


inline double VxcRef(const double val) {
	const double a = -0.61, b = 0.39, c = 0.33;
	const double v = abs2ToFilling(val);
	return a*(1+b)*pow(v,b) - c*v;
}

inline double ExcRef(const double val) {
	const double a = -0.61, b = 0.39, c = 0.33;
	const double v = abs2ToFilling(val);
	return (a*pow(v,b) - c*v/2)*val;
}


inline double VxcPaper(const double val) {
	const double a = -0.782133, b = 0.5, c = 0.33, f = 0.2774;
	const double v = abs2ToFilling(val);
	return a*(1+b)*pow(v,b) - c*v + 2*f*v;
}

inline double ExcPaper(const double val) {
	const double a = -0.782133, b = 0.5, c = 0.33, f = 0.2774;
	const double v = abs2ToFilling(val);
	return (a*pow(v,b) - c*v/2 + f*v)*val;
}


inline double VxcPaperRealFilling(const double v) {
	const double a = -0.782133, b = 0.5, c = 0.33, f = 0.2774;
	return a * (1 + b) * pow(v, b) - c * v + 2 * f * v;
}

inline double ExcPaperRealFilling(const double val, const double deltaB) {
	const double a = -0.782133, b = 0.5, c = 0.33, f = 0.2774;
	const double v = abs2ToFillingRealB(val, deltaB);
	return (a * pow(v, b) - c * v / 2 + f * v) * val;
}

/* Add XC potential with gradient dependence. This one is good because it is derived from the full energy of XC.
It is different from Vxcgrad in that it generates vp as a variable and use the gradient only when gradcoe != 0.0
I did not write a new and matching Exc with gradient bacause I actually do not use gradient correction now. */
void CalcVxcGradRealFilling(arma::mat& vp, const arma::mat& val) {
	iint ix, iy, jx, jy, kx, ky;
	const auto ax2 = ax * ax, ay2 = ay * ay;
	for (iint i = 0; i < Nxy; i++) {
		gridFromIndex(i, ix, iy);
		vp(ix, iy) = VxcRealFilling(val(ix, iy));
		const auto l = lnn(i), r = rnn(i), u = unn(i), d = dnn(i); // indices of neighbors. ==-1 if not found.

		/* find gradient of density from neighbours */
		if (gradcoe != 0.0) {
			if (l >= 0 && r >= 0) {
				gridFromIndex(l, jx, jy);
				gridFromIndex(r, kx, ky);
				vp(ix, iy) -= gradcoe * 0.25 / ax2 * pow(val(jx, jy) - val(kx, ky), 2);
			}
			else if (l >= 0) {
				gridFromIndex(l, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ax2 * pow(val(jx, jy) - val(ix, iy), 2);
			}
			else if (r >= 0) {
				gridFromIndex(r, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ax2 * pow(val(jx, jy) - val(ix, iy), 2);
			}
			if (u >= 0 && d >= 0) {
				gridFromIndex(u, jx, jy);
				gridFromIndex(d, kx, ky);
				vp(ix, iy) -= gradcoe * 0.25 / ay2 * pow(val(jx, jy) - val(kx, ky), 2);
			}
			else if (u >= 0) {
				gridFromIndex(u, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ay2 * pow(val(jx, jy) - val(ix, iy), 2);
			}
			else if (d >= 0) {
				gridFromIndex(d, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ay2 * pow(val(jx, jy) - val(ix, iy), 2);
			}
			/*laplacian*/
			vp(ix, iy) += gradcoe * 2.0 * pow(val(ix, iy), 2) * (2.0 / ax2 + 2.0 / ay2);
			if (l >= 0) {
				gridFromIndex(l, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * val(ix, iy) * val(jx, jy) / ax2;
			}
			if (r >= 0) {
				gridFromIndex(r, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * val(ix, iy) * val(jx, jy) / ax2;
			}
			if (u >= 0) {
				gridFromIndex(u, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * val(ix, iy) * val(jx, jy) / ay2;
			}
			if (d >= 0) {
				gridFromIndex(d, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * val(ix, iy) * val(jx, jy) / ay2;
			}
		}
		vp(ix, iy) *= scaleXC;//rescale the Vxc
	}
}

/* Add XC potential with gradient dependence. This one is not good because it is not derived from the full energy of XC.
It is different from Vxcgrad in that it generates vp as a variable*/
void CalcVxc(arma::mat& vp, const arma::mat& val) {
	//vp.zeros;//set all elements to zero
	iint ix, iy, jx, jy, kx, ky;
	const auto ax2 = ax * ax, ay2 = ay * ay;
	for (iint i = 0; i < Nxy; i++) {
		gridFromIndex(i, ix, iy);
		vp(ix, iy) = Vxc(val(ix, iy));// overwrite old values by the local part of Vxc

		const auto l = lnn(i), r = rnn(i), u = unn(i), d = dnn(i); // indices of neighbors. ==-1 if not found.
		/* find gradient of density from neighbours */
		if (gradcoe!=0.0){
		if (l >= 0 && r >= 0) {
			gridFromIndex(l, jx, jy);
			gridFromIndex(r, kx, ky);
			vp(ix, iy) += gradcoe * 0.25 / ax2 * pow(abs2ToFilling(val(jx, jy) - val(kx, ky)), 2);
		}
		else if (l >= 0) {
			gridFromIndex(l, jx, jy);
			vp(ix, iy) += gradcoe * 1.0 / ax2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
		}
		else if (r >= 0) {
			gridFromIndex(r, jx, jy);
			vp(ix, iy) += gradcoe * 1.0 / ax2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
		}
		if (u >= 0 && d >= 0) {
			gridFromIndex(u, jx, jy);
			gridFromIndex(d, kx, ky);
			vp(ix, iy) += gradcoe * 0.25 / ay2 * pow(abs2ToFilling(val(jx, jy) - val(kx, ky)), 2);
		}
		else if (u >= 0) {
			gridFromIndex(u, jx, jy);
			vp(ix, iy) += gradcoe * 1.0 / ay2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
		}
		else if (d >= 0) {
			gridFromIndex(d, jx, jy);
			vp(ix, iy) += gradcoe * 1.0 / ay2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
		}
		}
		vp(ix, iy) *= scaleXC;//rescale the Vxc
	}
}


/* Add XC potential with gradient dependence. This one is good because it is derived from the full energy of XC.
It is different from Vxcgrad in that it generates vp as a variable and use the gradient only when gradcoe != 0.0
I did not write a new and matching Exc with gradient bacause I actually do not use gradient correction now. */
void CalcVxcG(arma::mat& vp, const arma::mat& val) {
	iint ix, iy, jx, jy, kx, ky;
	const auto ax2 = ax * ax, ay2 = ay * ay;
	for (iint i = 0; i < Nxy; i++) {
		gridFromIndex(i, ix, iy);
		vp(ix, iy) = Vxc(val(ix, iy));

		const auto l = lnn(i), r = rnn(i), u = unn(i), d = dnn(i); // indices of neighbors. ==-1 if not found.
		/* find gradient of density from neighbours */
		if (gradcoe != 0.0) {
			if (l >= 0 && r >= 0) {
				gridFromIndex(l, jx, jy);
				gridFromIndex(r, kx, ky);
				vp(ix, iy) -= gradcoe * 0.25 / ax2 * pow(abs2ToFilling(val(jx, jy) - val(kx, ky)), 2);
			}
			else if (l >= 0) {
				gridFromIndex(l, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ax2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
			}
			else if (r >= 0) {
				gridFromIndex(r, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ax2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
			}
			if (u >= 0 && d >= 0) {
				gridFromIndex(u, jx, jy);
				gridFromIndex(d, kx, ky);
				vp(ix, iy) -= gradcoe * 0.25 / ay2 * pow(abs2ToFilling(val(jx, jy) - val(kx, ky)), 2);
			}
			else if (u >= 0) {
				gridFromIndex(u, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ay2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
			}
			else if (d >= 0) {
				gridFromIndex(d, jx, jy);
				vp(ix, iy) -= gradcoe * 1.0 / ay2 * pow(abs2ToFilling(val(jx, jy) - val(ix, iy)), 2);
			}
			/*laplacian*/
			vp(ix, iy) += gradcoe * 2.0 * pow(abs2ToFilling(val(ix, iy)), 2) * (2.0 / ax2 + 2.0 / ay2);
			if (l >= 0) {
				gridFromIndex(l, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * abs2ToFilling(val(ix, iy)) * abs2ToFilling(val(jx, jy)) / ax2;
			}
			if (r >= 0) {
				gridFromIndex(r, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * abs2ToFilling(val(ix, iy)) * abs2ToFilling(val(jx, jy)) / ax2;
			}
			if (u >= 0) {
				gridFromIndex(u, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * abs2ToFilling(val(ix, iy)) * abs2ToFilling(val(jx, jy)) / ay2;
			}
			if (d >= 0) {
				gridFromIndex(d, jx, jy);
				vp(ix, iy) -= gradcoe * 2.0 * abs2ToFilling(val(ix, iy)) * abs2ToFilling(val(jx, jy)) / ay2;
			}
		}
		vp(ix, iy) *= scaleXC;//rescale the Vxc
	}
}

/* Add XC potential. */
template<typename T>
void addXC(T& ham, const arma::mat &val) {
	iint ix, iy;
	for(iint i=0; i<Nxy; i++) {
		gridFromIndex(i, ix, iy);
		ham[At(i,i)] += Vxc(val(ix,iy))*scaleXC;
	}
}

#endif
