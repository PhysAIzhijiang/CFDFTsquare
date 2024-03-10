#!/bin/sh
ARMA_DIR=/sharedata01/yhu/workspace/data/CFDFTsquare/SymFractionTrans/arma
LD_LIBRARY_PATH=$ARMA_DIR/lib64:$LD_LIBRARY_PATH $*
