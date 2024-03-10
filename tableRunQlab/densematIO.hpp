#ifndef DENSEMATIO_HPP_
#define DENSEMATIO_HPP_

#include <cstdint>
#include <complex>
#include <string>
#include <iostream>
#include <armadillo>//error if deleted

/* 
 * (I)O for complex dense matrices.
 * Format:
 * Int64: version no.
 * type in 1,2,3 :: Int,  [Int,Float64,ComplexF64]
 * n_dims :: Int
 * Int64*2: nrows, ncols
 * ComplexF64(*): 
 */

std::int64_t writenumarray(std::string file, const arma::cx_mat mm)
{
    using namespace std;
    const int64_t ver = 1, kind = 3;
    int64_t i, j, ndim, ret;
    ofstream fh(file, ios::binary);
    i = mm.n_rows;
    j = mm.n_cols;
    ndim = 2;
    fh.write((char *)&ver, sizeof(ver));
    fh.write((char *)&kind, sizeof(kind));
    fh.write((char *)&ndim, sizeof(ndim));
    fh.write((char *)&i, sizeof(i));
    fh.write((char *)&j, sizeof(j));
    auto ptr = mm.memptr();
    ret = sizeof(int64_t) * 5;
    if (mm.n_elem > 0)
    {
        size_t len = sizeof(ptr[0]) * mm.n_elem;
        fh.write((char *)ptr, len);
        ret += len;
    }
    fh.close();
    return ret;
}

int readnumarray(std::string file, arma::cx_mat &mm){
    using namespace std;
    int64_t ver , kind;
    int64_t i, j, ndim, ret;
    ifstream fh(file, ios::binary);

    
    fh.read((char *)&ver, sizeof(ver));
    fh.read((char *)&kind, sizeof(kind));
    fh.read((char *)&ndim, sizeof(ndim));
    cout << ver << ' ' << kind << ' ' << ndim << endl;
    if (ver != 1 || kind != 3 || ndim != 2)
    {
        return -1;
    }

    fh.read((char *)&i, sizeof(i));
    fh.read((char *)&j, sizeof(j));
    mm.set_size(i, j);
    auto ptr = mm.memptr();
    std::cout << mm.n_elem << std::endl;
    ret = sizeof(int64_t) * 5;
    if (mm.n_elem > 0)
    {
        size_t len = sizeof(ptr[0]) * mm.n_elem;
        fh.read((char *)ptr, len);
        ret += len;
    }
    fh.close();
    return ret;
}

#endif
