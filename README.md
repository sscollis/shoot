# Shoot

Compressible stability solver using shooting.

**NOTE:** this has now been tested against `stab` and gives identical 
polished results.

## Build

```bash
ln -s gcc.mak Makefile
make USE_NR=1
```
Note that there are three options in building depending on the type of
mean flow file that you wish to use (this should really be a runtime
option).  Set `MEAN_2D=1` to use an LNS file for the mean flow profile
(this should be on a body-fitted mesh) and use `MEAN_BL=1` to use a BL
file.  The default is to read the mean flow from a collection of 
profile (note that `npost` can generate profiles from an LNS file).

## Running

That said, `shoot` is quite flexible and allow you to polish eigenvalues, 
solve the adjoint, and compute nonparallel terms all for a variety of 
mean flows and formats.

The `thesis` test case in `stab` exercises `shoot` to polish the spatial
eigensolution and output the regular and adjoint eigenfuncations.

The idea is that `shoot` complements `stab` by allowing you to polish, compute
adjoints, and include nonparallel effects.

## Notes
1. Currently this uses the Numerical Recipes RTSAFE routine (not included)
   so that you need to provide that (or implement another root finder)
2. I like `zeroin` function that is publically available on Netlib and
   that would be trivial to implement.
3. I have only tested the forward and adjoint solvers, not yet the nonparallel
   terms so please USE WITH CAUTION.

S. Scott Collis\
flow.physics.simulation@gmail.com
