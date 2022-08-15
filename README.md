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
mean flow file that you wish to use and there are three different versions of 
the code built for each type of mean flow: `shoot` works with a collection of 
2d profiles.  `shoot-2d` uses a special body-fitted mesh generated using 
`npost -l` when run on a LNS3D mean solution file. Finally, there is an experimental
version called `shoot-bl` that uses a BL format from NASA. 
The standard quasi-parallel linear stability analysis generally uses `shoot` while 
You must use `shoot-2d` to incorporate nonparallel effects.  NOTE: that `npost -p` can generate profiles from an LNS mean flow file for use by `shoot`. 

## Running

That said, `shoot` is quite flexible and allow you to polish eigenvalues, 
solve the adjoint, and compute nonparallel terms all for a variety of 
mean flows and formats.

The `thesis` test case in `stab` exercises `shoot` to polish the spatial
eigensolution and output the regular and adjoint eigenfuncations.

The idea is that `shoot` complements `stab` by allowing you to polish, compute
adjoints, and include nonparallel effects.

That case is in the `test` directory and can be run using:
```bash
cd test
./run.sh
```
Notes:
  1. This uses `ndiff` which must be in your path to do a numerical
     difference of the regular and adjoint eigenfunctions to make sure that
     there is no regression.
  2. To run without the regression test, enter `../shoot.exe < shoot.exe`
  3. You can visualize the regular and adjoint eigenfunctions using Gnuplot 
     with the `efun.com` and `adj.com` scripts.

## Notes
1. Currently this uses the Numerical Recipes RTSAFE routine (not included)
   so that you need to provide that (or implement another root finder)
2. I like `zeroin` function that is publically available on Netlib and
   that would be trivial to implement.
3. I have tested the forward, adjoint solvers and nonparallel correction
   terms.  See the `lns3d/test/pcyl` case for an example. 
   
## Tests 
 
The test runs for `shoot` are found in the `stab` and `lns3d` repositories. 

S. Scott Collis\
flow.physics.simulation@gmail.com
