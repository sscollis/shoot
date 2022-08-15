# Tollmein-Schilichting Test case 

This is the TS test case from Collis Ph.D. Thesis, Chapter 4.
It is computing using the `fsc` compressible Falkner-Skan solver 
for the mean flow which matches the thesis results exactly.

## Spatial results

```bash
FSC: M=0.3, lambda=0, beta_h=0, Tw/T0=1, Re_\delta_1 = 1000, Pr=1

The polished eigenvalue is:

alpha = k_x = 2.2804739839207E-01 -6.5163197601963E-03

\omega = 0.08
\lambda_{TS} = 2\pi/k_x = 27.552103068969444
```

## Running the test

The script `run.sh` polishes from the given eigenvalue and outputs the regular
and adjoint eigenfunction.
```bash
./run.sh
```
## Plotting the results

Plot the regular eigenfunction (normalized) using
```bash
gnuplot
load 'eig.com'
```
and similarly the adjoint eigenfunction using `adj.com`.  

## Notes
1. The orthogonality between regular and adjoint is not perfect as the output
   is only an approximation using trapezoidal integration whereas we should 
   use RK4 to be consistent.
2. 


S. Scott Collis\
flow.physics.simulation@gmail.com
