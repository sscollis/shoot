#!/bin/bash
#
# Make for nonparallel
#
#make clean && make USE_NR=1 DEFINES=-DUSE_BSLIB MEAN_2D=1
#
# Make for parallel
#
#make clean && make USE_NR=1 DEFINES=-DUSE_BSLIB
#
# this one works
#
make clean && make USE_NR=1
exit 0
