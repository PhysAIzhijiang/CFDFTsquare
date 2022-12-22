// This file is under public domain.
#ifndef COMMONINCLUDE_HPP
#define COMMONINCLUDE_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <type_traits>
#include <chrono>
#include <string>

#include <utility>
#include <map>
#include <unordered_map>

#include <array>
#include "tools.hpp"

#include<algorithm> // for copy
#include <iomanip> // for std::setprecision, writing data into files

typedef int iint; // armadillo want int, for matrix multiplication (BLAS), etc. So does mkl_sparse
typedef std::complex<double> ComplexF64;
/* custom universe */
const double Pi = acos(-1);
const ComplexF64 I {0,1};

typedef std::pair<iint,iint> At;


#endif
