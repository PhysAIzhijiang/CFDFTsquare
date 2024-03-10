// This file is under public domain.
#ifndef EIGENSOLVER_HPP
#define EIGENSOLVER_HPP
//header guards to prevent multiple inclusions of contents.
/*
 * Various sparse matrix eigensolvers
 */

#include <iostream>
#include <cstdio>
#include <cstring>
/*
 * MKL FEAST
 */
#include <mkl.h>

/* feastinit (MKL_INT* fpm);
 * [0]: verbose?
 * [4]: initial subspace?
 */

namespace FEAST {
	MKL_INT fpm[128], loop;
	double epsout, *res = nullptr; //C++ feature: declare two variables of the same type on the same line;
	MKL_INT lres = 0;// length of array res, should be kept track of to ensure lres=m0 on output
	MKL_INT info = 0;
};

void mkl_feast_init() {
	feastinit(FEAST::fpm);
}

void mkl_feast_prepare(MKL_INT Norb) {
	if (FEAST::lres < Norb) {
		FEAST::lres = Norb;// lres keeps track length of array res
		if (FEAST::res != nullptr) {// What about the zeroth run? res is nullptr, is that OK?
			delete[] FEAST::res;
		}
		FEAST::res = new double[FEAST::lres];//???????
	}
}

/*
 * void zfeast_hcsrev (const char * uplo, const MKL_INT * n, const MKL_Complex16 * a,
 const MKL_INT * ia, const MKL_INT * ja, MKL_INT * fpm, double * epsout, MKL_INT * loop,
 const double * emin, const double * emax, MKL_INT * m0, double * e, MKL_Complex16 * x,
 MKL_INT * m, double * res, MKL_INT * info);
 *
 * Warning: Entries in ia, ja must be >=1.
 *
 * Info for info:
 * 202	Error	Problem with size of the system n (n≤0)
 * 201	Error	Problem with size of initial subspace m0 (m0≤0 or m0>n)
 * 200	Error	Problem with emin,emax (emin≥emax)
 * (100+i)	Error	Problem with i-th value of the input Extended Eigensolver parameter (fpm[i - 1]). Only the parameters in use are checked.
 * 4	Warning	Successful return of only the computed subspace after call with
 * fpm[13] = 1
 * 3	Warning	Size of the subspace m0 is too small (m0<m)
 * 2	Warning	No Convergence (number of iteration loops > fpm[3])
 * 1	Warning	No eigenvalue found in the search interval. See remark below for further details.
 * 0	Successful exit
 * -1	Error	Internal error for allocation memory.
 * -2	Error	Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input.
 * -3	Error	Internal error of the reduced eigenvalue solver Possible cause: matrix B may not be positive definite. It can be checked by setting fpm[27] = 1 before calling an Extended Eigensolver routine, or by using LAPACK routines.
 * -4	Error	Matrix B is not positive definite.
 */

const char *mkl_eigen_errors[] = {
	"***	Error	Unknown error.",
	"(100+i)	Error	Problem with i-th value of the input Extended Eigensolver parameter (fpm[i-1]). Only the parameters in use are checked.",
	"202	Error	Problem with size of the system n (n≤0)",
	"201	Error	Problem with size of initial subspace m0 (m0≤0 or m0>n)",
	"200	Error	Problem with emin,emax (emin≥emax)",
	"(100+i)	Error	Problem with i-th value of the input Extended Eigensolver parameter (fpm[i-1]). Only the parameters in use are checked.",
	"4	Warning	Successful return of only the computed subspace after call with fpm[13] = 1",
	"3	Warning	Size of the subspace m0 is too small (m0<m)",
	"2	Warning	No Convergence (number of iteration loops > fpm[3])",
	"1	Warning	No eigenvalue found in the search interval. See remark below for further details.",
	"0	Successful exit",
	"-1	Error	Internal error for allocation memory.",
	"-2	Error	Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input.",
	"-3	Error	Internal error of the reduced eigenvalue solver Possible cause: matrix B may not be positive definite. It can be checked by setting fpm[27] = 1 before calling an Extended Eigensolver routine, or by using LAPACK routines.",
	"-4	Error	Matrix B is not positive definite.",
};

const char *feast_error_info(const MKL_INT ierr) {
	using namespace std;
	char serr[13];

	snprintf(serr, 13, sizeof(ierr)==4 ? "%d" : "%ld", ierr);
	const int len = strlen(serr);

	/*(100+i) errors */
	if(ierr>=100 && ierr<200) {
		return mkl_eigen_errors[1];
	}

	/* simple errors */
	for(auto s : mkl_eigen_errors) {
		if(! strncmp(serr,s,len)) {
			return s;
		}
	}

	return mkl_eigen_errors[0];
}

MKL_INT mkl_feast_ev(const char uplo, const MKL_INT n,
		const MKL_Complex16 * a, const MKL_INT * ia, const MKL_INT * ja,
		double emin, double emax, MKL_INT Norb, double *ees, MKL_Complex16 *evs,
		bool useGuess) {

	MKL_INT nee;/* This is the return value, the number of eigenstates found in the energy interval*/
	mkl_feast_prepare(Norb);
	FEAST::fpm[4] = useGuess ? 1 : 0;
	if(*ia <= 0) {
		std::cerr << "mkl_feast_ev: Entries in ia, ja must be >=1." << std::endl;
		exit(-1);
	}

	zfeast_hcsrev(&uplo, &n, a, ia, ja, FEAST::fpm,//Array, dimension of 128. This array is used to pass various parameters to Extended Eigensolver routines.
		&FEAST::epsout,/*On output, contains the relative error on the trace: |tracei - tracei-1| /max(| emin|, |emax|)*/
		&FEAST::loop, /*On output, contains the number of refinement loop executed. Ignored on input.*/
		&emin,&emax,
		&Norb, /* On entry, specifies the initial guess for subspace dimension to be used, 0<m0<=n. Set m0>=m where m is the total number of eigenvalues located in the interval [emin, emax]. If the initial guess is wrong, Extended Eigensolver routines return info=3. */
		ees, /* Array of length m0. On output, the first m entries of e are eigenvalues found in the interval.*/
		evs,/* On entry, if fpm[4]=1, the array x of size n by m contains a basis of guess subspace where n is the order of the input matrix.*/
			/*On output, the first m columns of x contain the orthonormal eigenvectors corresponding to the computed eigenvalues e, with the i-th column of x holding the eigenvector associated with e[i].*/
		&nee,/* The total number of eigenvalues found in the interval [emin, emax]: 0<=m<=m0 */
		FEAST::res,/* Array of length m0. On exit, the first m components contain the relative residual vector */
		&FEAST::info);/*If info=0, the execution is successful.*/
	if(FEAST::info != 0) {
		std::cerr << "MKL FEAST error: " << FEAST::info << std::endl;
	}
	return nee;
}

void mkl_skip_feast() {
	FEAST::fpm[13] = 1;/*0:Standard use for Extended Eigensolver routines;1:Non-standard use for Extended Eigensolver routines: return the computed eigenvectors subspace after one single contour integration.*/
}

void mkl_feast_verbose() {
	FEAST::fpm[0] = 1;/*Specifies whether Extended Eigensolver routines print runtime status.0:Extended Eigensolver routines do not generate runtime messages at all;1:Extended Eigensolver routines print runtime status to the
	screen.*/
}
#endif
