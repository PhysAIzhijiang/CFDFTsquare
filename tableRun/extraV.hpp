// This file is under public domain.
#ifndef EXTRAV_HPP
#define EXTRAV_HPP

#include "commoninclude.hpp"
#include "sparsemat.hpp"
#include "lattice.h"
#include "syspara.hpp"
#include <armadillo>
#include <algorithm>    // std::max

void uniformdis(arma::mat& dis) {
	arma::arma_rng::set_seed_random();  // set the seed to a random value
	dis = dis.randu() - 0.5;
}

/*long-range impurity, impurity density or impurity number, strength, range*/
void longRangeImpurity(arma::mat& dis) {
	dis.zeros();
	const double mx = (Nx - 1) / 2.0, my = (Ny - 1) / 2.0;
	const double h = 1.0;
	for (iint i = 0; i < Nx; i++) {
		for (iint j = 0; j < Ny; j++) {
			double r2 = (i - mx) * (i - mx) * ax * ax * DielecX_vs_DielecY + (j - my) * (j - my) * ay * ay + pow(h,2);
			r2 = sqrt(r2);
			dis(i, j) = exp(-r2/2.0) / r2;
		}
	}
}

/* Set a Rectangular region connecting from the center to the edge to be a bump potential. */
void linedis(arma::mat& dis) {
	dis.zeros();
	const int mx = round((Nx - 1) / 2.0), my = round((Ny - 1) / 2.0);
	const int xwidth = Nrx;
	for (iint i = mx - xwidth; i < mx + xwidth; i++) {
		for (iint j = my; j < Ny; j++) {
			dis(i, j) = 1.0;
		}
	}
}


/* Set a Rectangular region connecting from the center to the edge to be a bump potential. */
void addDeltaBarrier(arma::mat& DeltaBarrierPotential, const enum SHAPE shapeCS, const double DeltaBarrierNp, const double geometryYoverXDeltaBarrier) {
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

}


/* Set a Rectangular region of density to `filling`. `sideX` is in number of grids. */
/*
void fillRectangle(arma::mat& density, const double filling, double sideX, double sideY) {
	sideX = sideX / 2.0;
	sideY = sideY / 2.0;
	const double wfsq = fillingToAbs2(filling);
	const double mx = (Nx - 1) / 2.0, my = (Ny - 1) / 2.0;
	for (iint j = 0; j < Ny; j++) {
		for (iint i = 0; i < Nx; i++) {
			const double di = i - mx, dj = j - my;
			if (abs(di) <= sideX && abs(dj) <= sideY) {
				density(i, j) = wfsq;
			}
			else {
				density(i, j) = 0;
			}
		}
	}
}
*/




template<typename iint, typename T>
void diskWell(MMSpmatR<iint,T> &ham, double boundHeight) {
	const double mx = (Nx-1) / 2.0, my = (Ny-1) / 2.0;
	const double r = min(mx,my);	
	for(iint i = 0; i<Nx; i++) {
		for(iint j=0; j<Ny; j++) {
			if((i-mx)*(i-mx) + (j-my)*(j-my) > r*r) {
				auto ix = gridToIndex(i,j);
				ham[At(ix,ix)] += boundHeight; 
			}
		}
	}
}

template<typename iint, typename T>
void noop(MMSpmatR<iint,T> &ham, double boundHeight) {
}


template<typename T>
void analyticdisk(T& Vdisk, double boundHeight) {
	const double mx = (Nx - 1) / 2.0, my = (Ny - 1) / 2.0;
	const double r = min(mx, my);
	for (iint i = 0; i < Nx; i++) {
		for (iint j = 0; j < Ny; j++) {
			//if ((i - mx) * (i - mx) + (j - my) * (j - my) > r * r) {
				//auto ix = gridToIndex(i, j);
				Vdisk(i,j) = boundHeight;
			//}
		}
	}
}


#endif
