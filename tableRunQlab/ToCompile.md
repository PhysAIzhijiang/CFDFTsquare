# CFDFT-square manual

**Author: Arvin-Hu**

## System environment


(to load intel mkl impi)

module load devtoolset/gcc-11 intel/2022 
make arma
make cfdft-unopt

(if gcc version not correct, there will be errors like)

/opt/rh/devtoolset-11/root/usr/include/c++/11/bits/stringfwd.h(72): error: expected a ";"
    template<typename _CharT, typename _Traits = char_traits<_CharT>, 

/opt/rh/devtoolset-11/root/usr/include/c++/11/bits/stl_iterator.h(1120): error: invalid combination of type specifiers
      inline bool


(to calculate)

vi syspara.hpp

mkdir 20x20_12

cd 20x20_12

cp ../cfdft.x .  put excutable here

./cfdft.x  make sure excutable can run

../qopen 20 16g ./cfdft.x submit job. Sometimes need chmod +x ../qopen


**About makefile**

1. 第一个MKLROOT在你加载Intel MKL库”在qlab集群上执行module load mkl“之后环境变量里面就有了，makefile会自动找到的，这边“？=”的意思是环境变量未设置时的值。当然你也可以在加载MKL库之后通过“echo $MKLROOT"来查看其值。

2. IMPIROOT 这个变量整个文件中没有使用，所以我不完全确定这边要的值是什么形式。但应该就是加载intel MPI模块之后”module load mpi"之后的环境变量I_MPI_ROOT "echo $I_MPI_ROOT"的值(/opt/intel/oneapi/mpi/2021.6.0)。

3. OpenMPIROOT这个变量由于CentOS和debian(这个makefile文件是针对debian系统写的)的文件位置不同，建议直接将“-I${OpenMPIROOT}/include” 替换为“-I/usr/include/openmpi-x86_64”

4. CPLUS_INCLUDE_PATH 一般是这个/usr/include/c++/4.8.5/，可能需要新版的gcc，那会在另外一个地方



**Qlab not installed**

PETSC_DIR和SLEPC_DIR我们集群上目前没有


