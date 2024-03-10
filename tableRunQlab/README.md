
# This file is under public domain.

## Dependences
Require Intel(r) MKL and intel C++ compiler to compile.
The program uses armadillo library, http://arma.sourceforge.net, for convenience.

## Files

commoninclude.hpp        Header files and basic variable types.

extraV.hpp	         Extra external potential *

DFTbuilder.hpp		 Hamiltonian builder *

eigensolver.hpp		 Matrix diagonalization, of eigenval in some interval

electromag.hpp		 A & V potential, using convolution

excorr.hpp		 Exchange correlation functions. *

lattice.h		 Lattice builder.

main.cpp		 Main program.

mklspmat.hpp		 Sparse matrix multiplication.

sparsemat.hpp		 Sparse matrix types.

syspara.hpp		 System parameters and settings. **

teestream.hpp		 Log file writing utilities.

main_util.hpp		 Record generation, output file writing, etc. Used by main in namespace MAIN. */2

tools.hpp		 Tools for c++ programming.

type.hpp		 Tools for c++ types.


## Tips
The armadillo library can be built & installed locally with `make arma`. Remember to remove the arma folder before doing so to ensure a fresh make.
Then `make cfdft` compiles the program, according to parameters in "syspara.hpp", etc.
