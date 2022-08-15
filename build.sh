#!/bin/bash
#
# Build a fresh version -- note that threee different versions are built
# depending on the type of base flow you intend to use:
#
# shoot.exe   :  collection of 2d profiles
# shoot-2d.exe:  LNS meanflow on body-fitted mesh (for nonparallel effects)
# shoot-bl.exe:  Streett's BL format (experimental -- not currently working)
#
# Revised:  8/15/22
# Author:   S.Scott Collis 
#
set -e
make clean && make USE_NR=1 $@
exit $? 
