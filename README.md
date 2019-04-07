# EZ-AFM
Quick and easy simulations of atomic force microscopy (AFM)-like pulling of macromolecules, using a simple Go-like model with Brownian dynamics.

## Potential of the system
The total configuration potential is the sum of the following energy terms:

![equation 1](demo/eqn/eqn1.png)

where:

![equation 2](demo/eqn/eqn2.png)

![equation 3](demo/eqn/eqn3.png)

![equation 4](demo/eqn/eqn4.png)

![equation 5](demo/eqn/eqn5.png)

## Brownian Dynamics without hydrodynamics interaction
The displacement propagation of the system is obtained by solving the following equation numerically:

![equation 6](demo/eqn/eqn6.png)

## Brownian Dynamics with hydrodynamics interaction 
#### (under testing and debugging)

When considering hydrodynamics interactions, the displacement propagation can be expressed as:

![equation 7](demo/eqn/eqn7.png)

where **D** is the Ronte-Prager-Yamakawa diffusion tensor:

![equation 8](demo/eqn/eqn8.png)

**g** is a Gaussian distributed array and **S** is obtained by Cholesky decomposition of **D** since:

![equation 9](demo/eqn/eqn9.png)


## Notes
- Generate paramter files using gen-top, then run dynamics with main-prog.
- All the files (`crd`, `psf`, `dcd`, `prm` and `rtf`) generated will be compatible with CHARMM package for further analysis.
- Program currently under testing and debugging. 
- 3 short demos of ubiquitin are provided (Note the AFM-like simulations are pulling at high speed for quick demonstration):

1. Short equilibration

<img src="demo/equil.gif" width="50%" height="50%" />

2. AFM-like pulling with N terminal fixed

<img src="demo/fixed.gif" width="50%" height="50%" />

3. AFM-like pulling with out fixing any atom

<img src="demo/free.gif" width="50%" height="50%" />
