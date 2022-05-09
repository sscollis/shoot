# Shoot

Compressible stability solver using shooting

## Build

    ln -s gcc.mak Makefile
    make USE_NR=1

## Notes
1. Currently this uses the Numerical Recipes RTSAFE routine (not included)
   so that you need to provide that (or implement another root finder)

S. Scott Collis\
Mon May  9 09:33:41 MDT 2022
