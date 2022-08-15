#!/bin/bash
set -e
#
# Shooting to cleanup spatial 
#
../shoot.exe < shoot.inp
#
# check quality of solution 
#
echo "============================================="
echo "ndiff of efun" && \
ndiff -s -abserr 1e-6 efun.out efun.ref | tee efun.log  && \
echo "=============================================" && \
echo "ndiff of adj" && \
ndiff -s -abserr 1e-6 adj.out adj.ref | tee adj.log  && \
status=$?
echo "============================================="
echo 'ndiff completed with status:' $status
exit $status 
