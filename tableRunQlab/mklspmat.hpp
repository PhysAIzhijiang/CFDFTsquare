// This file is under public domain.
#ifndef MKLSPMAT_HPP
#define MKLSPMAT_HPP

#include "commoninclude.hpp"
#include <mkl_spblas.h>
#include "sparsemat.hpp"

void mkl_spblas_call_and_check(sparse_status_t rv) {
	using namespace std;
	if(rv != SPARSE_STATUS_SUCCESS) {
		std::cerr << "mkl sparse error code: " << rv << '.' << std::endl;
		std::cerr << "Aborting." << std::endl;
		exit(rv);
	}
}


template<typename iint>
struct MMSpmatR_op {
	MMSpmatR<iint,std::complex<double> > *spmat;
	struct matrix_descr descr;
	sparse_matrix_t csr;
	bool initialized = false;
};


template<typename iint>
struct MMSpmatR_op<iint> create_Hermitian_op(MMSpmatR<iint,std::complex<double> > &ham, int N_vc, int N_mm) {
	//add information that ham is hermitian to speed up
	auto cc = mkl_spblas_call_and_check;
	struct MMSpmatR_op<iint> matop;

	matop.spmat = &ham;
	cc( mkl_sparse_z_create_csr(&matop.csr, SPARSE_INDEX_BASE_ZERO,
				    matop.spmat->n_rows, matop.spmat->n_cols,
				    matop.spmat->rowstart, matop.spmat->rowstart+1,
				    matop.spmat->cols,
				    (MKL_Complex16 *)matop.spmat->vals) );
	matop.descr.type = SPARSE_MATRIX_TYPE_HERMITIAN;
	matop.descr.mode = SPARSE_FILL_MODE_UPPER;
	matop.descr.diag = SPARSE_DIAG_NON_UNIT;	
	
//	cc( mkl_sparse_set_mm_hint(matop.csr,SPARSE_OPERATION_NON_TRANSPOSE,matop.descr,SPARSE_LAYOUT_ROW_MAJOR,N_vc,N_mm+1) );

//	cc( mkl_sparse_optimize(matop.csr) );

	{
		/* Premise check:
		 * Internal structures of MKL sparse_matrix_t can be used freely. 
		 * Checking here. 
		 */
		auto rs = matop.spmat->rowstart, re = matop.spmat->rowstart+1;
		auto cols = matop.spmat->cols;
		auto vals = matop.spmat->vals;
		auto nr = matop.spmat->n_rows, nc = matop.spmat->n_cols;
		sparse_index_base_t idx;
		rs = nullptr; re=nullptr; cols = nullptr; vals = nullptr;
		
		cc( mkl_sparse_z_export_csr(matop.csr,&idx,&nr,&nc,&rs,&re,&cols,(MKL_Complex16**)&vals) );

		assert(cols == matop.spmat->cols);
		assert(rs == matop.spmat->rowstart);
		assert(re == matop.spmat->rowstart+1);
		assert(vals == matop.spmat->vals);
	}

	matop.initialized = true;
	return matop;
}


//put matop*x(:,1:xCol) into y
template<typename iint>
void mkl_mm(struct MMSpmatR_op<iint> matop, iint xCol, const MKL_Complex16 *x, MKL_Complex16 *y){
	auto cc = mkl_spblas_call_and_check;
	const std::complex<double> emc = 1.0, omc = 0.0;
	const MKL_Complex16 *emp = reinterpret_cast<const MKL_Complex16*>(&emc);
	const MKL_Complex16 *omp = reinterpret_cast<const MKL_Complex16*>(&omc);

	auto ld = matop.spmat->n_rows;

	if(! matop.initialized) {
		std::cerr << "matop uninitialized. Aborting." << std::endl;
		exit(1);
	}
	
	cc( mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE, *emp, matop.csr, matop.descr, SPARSE_LAYOUT_COLUMN_MAJOR, (const MKL_Complex16*) x, xCol, ld, *omp, (MKL_Complex16*) y, ld) );

	{
		/* Premise check:
		 * Internal structures of MKL sparse_matrix_t can be used freely. 
		 * Checking here. 
		 */
		auto rs = matop.spmat->rowstart, re = matop.spmat->rowstart+1;
		auto cols = matop.spmat->cols;
		auto vals = matop.spmat->vals;
		auto nr = matop.spmat->n_rows, nc = matop.spmat->n_cols;
		sparse_index_base_t idx;
		rs = nullptr; re=nullptr; cols = nullptr; vals = nullptr;
		
		cc( mkl_sparse_z_export_csr(matop.csr,&idx,&nr,&nc,&rs,&re,&cols,(MKL_Complex16**)&vals) );

		assert(cols == matop.spmat->cols);
		assert(rs == matop.spmat->rowstart);
		assert(re == matop.spmat->rowstart+1);
		assert(vals == matop.spmat->vals);
	}


}
#endif
