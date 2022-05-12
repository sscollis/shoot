# Shoot

Compressible stability solver using shooting.

**NOTE:** this has not been tested recently so please use with extreme
caution.

## Build

    ln -s gcc.mak Makefile
    make USE_NR=1

## Running

I really need to setup a test case for this.  That said, `shoot` is quite 
flexible and allow you to polish eigenvalues, solve the adjoint, and compute
nonparallel terms all for a variety of mean flows and formats.

## Notes
1. Currently this uses the Numerical Recipes RTSAFE routine (not included)
   so that you need to provide that (or implement another root finder)
2. I like `zeroin` function that is publically available on Netlib and
   that would be trivial to implement.
3. I have not tested this solver recently so please USE WITH CAUTION

--
S. Scott Collis
