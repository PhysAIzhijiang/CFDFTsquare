/* This file is under public domain. */
#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "commoninclude.hpp"

/*
 * Square spacial grid.
 */
iint unn(const iint);
iint lnn(const iint);
iint rnn(const iint);
iint dnn(const iint);
iint gridToIndex(const iint, const iint);
void gridFromIndex(const iint, iint&, iint&);

/* As noted before, 0-based index. */
const int BaseIndex = 0;
inline iint gridToIndex(const iint cx, const iint cy) {
	return (cy)*Nx + cx;
}
inline void gridFromIndex(const iint index, iint &cx, iint &cy) {
	cx = (index - BaseIndex) % Nx + BaseIndex;
	cy = (index - BaseIndex) / Nx + BaseIndex;
}


/* Indices to nearest neighbours. Return <0 if none. */
iint dnn(const iint idx) {
	iint m, n;
	gridFromIndex(idx, m, n);
	if(n>BaseIndex) {
		return gridToIndex(m,n-1);
	} else {
		return -1;
	}
}
iint lnn(const iint idx) {
	iint m, n;
	gridFromIndex(idx, m, n);
	if(m>BaseIndex) {
		return gridToIndex(m-1,n);
	} else {
		return -1;
	}
}
iint rnn(const iint idx) {
	iint m, n;
	gridFromIndex(idx, m, n);
	if(m+1 < Nx+BaseIndex) {
		return gridToIndex(m+1,n);
	} else {
		return -1;
	}
}
iint unn(const iint idx) {
	iint m, n;
	gridFromIndex(idx, m, n);
	if(n+1 < Ny+BaseIndex) {
		return gridToIndex(m,n+1);
	} else {
		return -1;
	}
}

#endif
