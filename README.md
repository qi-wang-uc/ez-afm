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

## Example input
```
!o CHARMM-like input but not exactly the same.
!o This program does not have the 4-character rule, so make sure use the full name of keywords.
!o Values in the triangle bracket are the default value.
!o Use curly brackets {} to refer a variable to eliminate ambiguity such as @a and @ab.

set protname = 1ubq	   ! The equal sign (=) is optional.
set nstep    = 10000      ! <100>
set tstep    = 1	   ! <1> timestep in ps.
set zeta     = 50	   ! <50>	damping coefficient.
set temp     = 298	   ! <298> temperature in Kelvin.
set outfreq  = 100 	   ! <10> frequency to write energy info and .afm files.
set dcdfreq  = 10	   ! <10> frequency to write coordinate
set nbdfreq  = 10	   ! <10> frequency to update nonbonded list
set dijfreq  = 10	   ! <10> frequency to update diffusion tensor, will be ignored when "hydro" is set to "none".
set hydro    = 0	   ! <0>  hydrodynamics (1) or none (0)

! Read toppar and coordinate into system
system prm    @{protname}.prm
system psf    @{protname}.psf
system cor    @{protname}.cor
! fix atom index
system fix    1

! AFM setup
afm nterm 1 cterm 76 force 2 velocity 0.000001 maxdist 280

! Run dynamics
dyna hydro @{hydro} nstep @{nstep} tstep @{tstep} temp @{temp} zeta @{zeta} -
     dcdfreq @{dcdfreq} nbdfreq @{nbdfreq} outfreq @{outfreq} dijfreq @{dijfreq} - 
     dcdname @{protname}.dcd

stop
```

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
