# EZ-AFM
Quick and easy simulations of atomic force microscopy (AFM)-like pulling of macromolecules, using a simple Go-like model with Brownian dynamics.

## Potential of the system
The total configuration potential is the sum of the following energy terms:

<a href="https://www.codecogs.com/eqnedit.php?latex=V_{\mathrm{configuration}}&space;=&space;V_{\mathrm{bonds}}&space;&plus;&space;V_{\mathrm{angles}}&space;&plus;&space;V_{\mathrm{dihedrals}}&space;&plus;&space;V_{\mathrm{non-bonded&space;}}" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_{\mathrm{configuration}}&space;=&space;V_{\mathrm{bonds}}&space;&plus;&space;V_{\mathrm{angles}}&space;&plus;&space;V_{\mathrm{dihedrals}}&space;&plus;&space;V_{\mathrm{non-bonded&space;}}" title="V_{\mathrm{configuration}} = V_{\mathrm{bonds}} + V_{\mathrm{angles}} + V_{\mathrm{dihedrals}} + V_{\mathrm{non-bonded }}" /></a>

where:

<a href="https://www.codecogs.com/eqnedit.php?latex=V_{\mathrm{bonds}}=\sum_{\mathrm{bonds}}&space;k_b(b-b_0)^2" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_{\mathrm{bonds}}=\sum_{\mathrm{bonds}}&space;k_b(b-b_0)^2" title="V_{\mathrm{bonds}}=\sum_{\mathrm{bonds}} k_b(b-b_0)^2" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=V_{\mathrm{angles}}=\sum_{\mathrm{angles}}&space;k_{\theta}(\theta-\theta_0)^2" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_{\mathrm{angles}}=\sum_{\mathrm{angles}}&space;k_{\theta}(\theta-\theta_0)^2" title="V_{\mathrm{angles}}=\sum_{\mathrm{angles}} k_{\theta}(\theta-\theta_0)^2" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=V_{\mathrm{dihedrals}}=\sum_{\mathrm{dihedrals}}&space;k_{\phi}[1&plus;\cos(n\phi-\delta)]" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_{\mathrm{dihedrals}}=\sum_{\mathrm{dihedrals}}&space;k_{\phi}[1&plus;\cos(n\phi-\delta)]" title="V_{\mathrm{dihedrals}}=\sum_{\mathrm{dihedrals}} k_{\phi}[1+\cos(n\phi-\delta)]" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=V_{\mathrm{non-bonded}}&space;=&space;\sum_{\mathrm{non-native}}\epsilon_1&space;\big(&space;\frac{\sigma_{i,j}}{r_{i,j}}\big)^{12}&space;&plus;&space;\sum_{\mathrm{native}}\epsilon_2\big[&space;\big(&space;\frac{\sigma_{i,j}}{r_{i,j}}\big)^{12}-\big(\frac{\sigma_{i,j}}{r_{i,j}}&space;\big)^6\big]" target="_blank"><img src="https://latex.codecogs.com/png.latex?V_{\mathrm{non-bonded}}&space;=&space;\sum_{\mathrm{non-native}}\epsilon_1&space;\big(&space;\frac{\sigma_{i,j}}{r_{i,j}}\big)^{12}&space;&plus;&space;\sum_{\mathrm{native}}\epsilon_2\big[&space;\big(&space;\frac{\sigma_{i,j}}{r_{i,j}}\big)^{12}-\big(\frac{\sigma_{i,j}}{r_{i,j}}&space;\big)^6\big]" title="V_{\mathrm{non-bonded}} = \sum_{\mathrm{non-native}}\epsilon_1 \big( \frac{\sigma_{i,j}}{r_{i,j}}\big)^{12} + \sum_{\mathrm{native}}\epsilon_2\big[ \big( \frac{\sigma_{i,j}}{r_{i,j}}\big)^{12}-\big(\frac{\sigma_{i,j}}{r_{i,j}} \big)^6\big]" /></a>

## Brownian Dynamics
The displacement propagation of the system is obtained by solving the following equation numerically:

<a href="https://www.codecogs.com/eqnedit.php?latex=\xi&space;\frac{\mathrm{d}r}{\mathrm{d}t}=-\frac{\partial&space;V_i(r_i)}{\partial&space;r_i}&plus;g_i(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\xi&space;\frac{\mathrm{d}r}{\mathrm{d}t}=-\frac{\partial&space;V_i(r_i)}{\partial&space;r_i}&plus;g_i(t)" title="\xi \frac{\mathrm{d}r}{\mathrm{d}t}=-\frac{\partial V_i(r_i)}{\partial r_i}+g_i(t)" /></a>

## Notes
- Generate paramter files using gen-top, then run dynamics with main-prog.
- All the files (`crd`, `psf`, `dcd`, `prm` and `rtf`) generated will be compatible with CHARMM package for further analysis.
- Program currently under testing and debugging. 
- 3 short demos of ubiquitin are provided:

1. Short equilibration
<img src="https://github.com/wangqi1990uc/ez-afm/blob/master/equil.gif" width="50%" height="50%" />

2. AFM-like pulling with N terminal fixed
<img src="https://github.com/wangqi1990uc/ez-afm/blob/master/fixed.gif" width="50%" height="50%" />

3. AFM-like pulling with out fixing any atom
<img src="https://github.com/wangqi1990uc/ez-afm/blob/master/free.gif" width="50%" height="50%" />


## TODO
- [ ] Better user variable handling
- [ ] Add hydrodynamic interaction tensor
