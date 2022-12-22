# This file is under public domain.

MKLROOT?=/opt/intel/mkl
IMPIROOT=/opt/intel/compilers_and_libraries/linux/mpi/intel64
OpenMPIROOT=/usr/lib/x86_64-linux-gnu/openmpi
ARMA_DIR?=${PWD}/arma
CPLUS_INCLUDE_PATH?=/usr/include/x86_64-linux-gnu/c++/8 # adhoc fix on some linux machines for icc problems
# slepc variables
PETSC_DIR=/usr/lib/petscdir/3.7.5/x86_64-linux-gnu-complex
SLEPC_DIR=/usr/lib/slepcdir/3.7.3/x86_64-linux-gnu-complex

# slepc compiler & linker options
SLEPC_LINKOPT= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
GCC_SLEPC_LINKOPT= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
SLEPC_CCOPT=-I${MKLROOT}/include -I${PETSC_DIR}/include -I${SLEPC_DIR}/include -I${OpenMPIROOT}/include


# mkl compiler and linker options
MKL_CCOPT=-I${MKLROOT}/include

INTEL_CCOPT=-I${ARMA_DIR}/include

MKL_LINKOPT=-L${ARMA_DIR}/lib -L${ARMA_DIR}/lib64 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

# slepc includes
MANSEC=EPS
#include ${SLEPC_DIR}/lib/slepc/conf/slepc_common


.PHONY: slepc gcc_slepc mkl tgz arma

slepc:
	rm -f $@.x
	icpc -o $@.x -O3 -g -debug full -std=c++14 ${SLEPC_CCOPT} -ftz -fp-model source -fp-model fast=1 -parallel -qopenmp main.cpp -larmadillo ${SLEPC_LINKOPT} -${CLINKER} ${SLEPC_EPS_LIB}
#	echo
#	echo 
#	echo g++
#	echo
#	echo
#	c++ -std=c++14 -Wall ${CCOPT} -m64  -fopenmp main.cpp -larmadillo ${LINKOPT} -${CLINKER} ${SLEPC_EPS_LIB}

gcc_slepc:
	rm -f $@.x
	mpic++ -o $@.x  -O3 -g -std=c++14 ${SLEPC_CCOPT} -m64  main.cpp -larmadillo ${SLEPC_LINKOPT} -${CLINKER} ${SLEPC_EPS_LIB}

arma_to_bin: arma_to_bin.cpp
	c++ -o $@ -I${ARMA_DIR}/include -L${ARMA_DIR}/lib -L${ARMA_DIR}/lib64 $< -larmadillo


cfdft:
	rm -f $@.x
	CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH} icpc  -o $@.x -O3 -g -debug full \
		-std=c++14 \
		${INTEL_CCOPT} \
		${MKL_CCOPT} \
		-ftz -fp-model source -fp-model fast=1 -parallel -qopenmp \
		main.cpp -larmadillo ${MKL_LINKOPT} \
		-Wl,-rpath,${ARMA_DIR}/lib64

cfdft-unopt:
	rm -f $@.x
	CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH} icpc  -o $@.x -O0 -g -debug full \
		-std=c++14 \
		${INTEL_CCOPT} \
		${MKL_CCOPT} \
		-ftz -fp-model source -fp-model fast=1 \
		main.cpp -larmadillo ${MKL_LINKOPT} \
		-Wl,-rpath,${ARMA_DIR}/lib64

arma:
	rm -fr arma arma*[0-9]
	echo ${ARMA_DIR}
	tar xJf arma*.tar.xz
	cd arma*[0-9] && CC=icc CXX=icpc cmake . -DDETECT_HDF5=false -DCMAKE_INSTALL_PREFIX:PATH=${ARMA_DIR} &&	CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH} make install
	echo '#!/bin/sh' > ld-arma.sh
	echo 'ARMA_DIR='${ARMA_DIR} >> ld-arma.sh
	echo 'LD_LIBRARY_PATH=$$ARMA_DIR/lib64:$$LD_LIBRARY_PATH' '$$*' >> ld-arma.sh
	chmod +x ld-arma.sh

tgz:
	tar --exclude-backups --exclude-vcs -czf cfdft-cc.tgz `ls *.cpp *.hpp *.h arma*.tar.xz *.txt readnumarray.m makefile 2>/dev/null`
