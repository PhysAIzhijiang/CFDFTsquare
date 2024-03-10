// This file is under public domain.
#ifndef SPARSEMAT_HPP
#define SPARSEMAT_HPP
/*
 * Sparse Matrix implementation.
 * Assuming dimensions and location of nonzero entries do not change. 
 */


#include <utility>
#include <map>
#include <unordered_map>
#include <cstdlib>
#include <vector>

/*
 * Sparse matrix, in matrix market and CSR format.
 * (As is CSR), prepared to column major form.
 * No destructor here.
 */
template<typename iint, typename T>
struct MMSpmatR {
	iint n_rows, n_cols;
	iint nnz; // number of filled entries (not necessary nonzero !!!)
	/* Everything has 0-based index, i.e. runs from 0 to nnz-1 */
	iint *rows;
	iint *cols;
	T *vals;
	iint *rowstart;
	/* 1-based index partners, for MKL, etc. */
	iint *rows1;
	iint *cols1;
	iint *rowstart1;
	std::map<std::pair<iint,iint>,iint> entryloc;
	//iint *diag;
	T& operator[](std::pair<iint,iint> idx) {
		return vals[ entryloc[idx] ];
	};
};


///*
// * MatMap: map wrapper for a matrix type T
// * T must have fields:
// * - nnz
// * - rows
// * - cols1
// * - vals
// */
//template<typename T>
//class MatMap {
//public:
//	T mat;
//	map<At,iint> entryloc;
//
//	MatMap() {
//	}
//	MatMap(T a) {
//		mat = a;
//		for(iint i = 0; i < a.nnz; i++) {
//			entryloc[At(mat.rows[i],mat.cols[i])] = i;
//		}
//	}
//	T& operator[](At idx) {
//		return mat.vals[entryloc[idx]];
//	}
//};
//

/*
  Collect map<array<iint,2>,T>, into a sparse matrix MMSpmatR.
  All entries are kept, including zeros.
*/

template<typename iint, typename T>
MMSpmatR<iint,T> mapToMMSpmatR(const iint m, const iint n, std::map<std::pair<iint,iint>,T> &draft) {
	using namespace std;
	
	/* 0-based index */
	MMSpmatR<iint,T> mat {m, n, static_cast<iint>(draft.size())};// armadillo features : .size() returns the total number of elements; int +- 2*10^9 should be large enough 
	iint i;

	mat.rows = new iint[mat.nnz];
	mat.cols = new iint[mat.nnz];
	mat.vals = new T[mat.nnz];
	mat.rowstart = new iint[m+1];

	// mat.diag = new iint[min(m,n)]{-1};

	// dump the map to temp arrays
	vector<iint> ti(mat.nnz), tj(mat.nnz);
	vector<T> tv(mat.nnz);

	i = 0;

	for(auto&& kv : draft) {
		ti[i] = (kv.first).first;
		tj[i] = (kv.first).second;
		tv[i] = kv.second;
		i++;
	}

	// order index
	iint *l = new iint[10*mat.nnz];
	for(iint i=0; i<mat.nnz;i++){
		l[i]=i;
	}

	// sorting indices
	// !! This is unnecessary for map. Keep for now.
	//sort(l, l+mat.nnz, [ti,tj](const iint a, const iint b){ return ti[a] < ti[b] || (ti[a]==ti[b] && tj[a]<tj[b]); });

	// put things in order
	// Also fill rowstart for CSR format ...
	// (and diag) nope
	iint nextr = 0;
	for (i = 0; i < mat.nnz; i++) {
		const auto ii = ti[l[i]], jj = tj[l[i]];
		const auto vv = tv[l[i]];
		mat.rows[i] = ii;
		mat.cols[i] = jj;
		mat.vals[i] = vv;
		mat.entryloc[At(ii,jj)] = i;
//		if(mat.rows[i] == mat.cols[i]) {
//			mat.diag[mat.rows[i]] = i;
//		}
		for (; nextr <= mat.rows[i]; nextr++) {
			mat.rowstart[nextr] = i;
		}
	}
	for (; nextr < m + 1; nextr++) {
		mat.rowstart[nextr] = i;
	}

	/* Assign 1-based index arrays. */
	mat.rows1 = new iint[mat.nnz];
	mat.cols1 = new iint[mat.nnz];
	mat.rowstart1 = new iint[m + 1];
	for (iint i = 0; i < mat.nnz; i++) {
		mat.rows1[i] = mat.rows[i] + 1;
		mat.cols1[i] = mat.cols[i] + 1;
	}
	for (iint i = 0; i < m + 1; i++) {
		mat.rowstart1[i] = mat.rowstart[i] + 1;
	}

	return mat;
}

#endif
