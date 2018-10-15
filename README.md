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

## Brownian Dynamics without hydrodynamics interaction
The displacement propagation of the system is obtained by solving the following equation numerically:

<a href="https://www.codecogs.com/eqnedit.php?latex=\xi&space;\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}t}=-\frac{\partial&space;V_i(r_i)}{\partial&space;r_i}&plus;g_i(t)" target="_blank"><img src="https://latex.codecogs.com/png.latex?\xi&space;\frac{\mathrm{d}r}{\mathrm{d}t}=-\frac{\partial&space;V_i(r_i)}{\partial&space;r_i}&plus;g_i(t)" title="\xi \frac{\mathrm{d}r}{\mathrm{d}t}=-\frac{\partial V_i(r_i)}{\partial r_i}+g_i(t)" /></a>

## Brownian Dynamics with hydrodynamics interaction 
#### (under testing and debugging)

When considering hydrodynamics interactions, the displacement propagation can be expressed as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{r_{t&plus;\Delta&space;t}}=\mathbf{r_{t}}&plus;\frac{\Delta&space;t}{k_BT}\mathbf{D\cdot&space;F}&plus;\sqrt{2\Delta&space;t}\mathbf{S\cdot&space;g}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\mathbf{r_{t&plus;\Delta&space;t}}=\mathbf{r_{t}}&plus;\frac{\Delta&space;t}{k_BT}\mathbf{D\cdot&space;F}&plus;\sqrt{2\Delta&space;t}\mathbf{S\cdot&space;g}" title="\mathbf{r_{t+\Delta t}}=\mathbf{r_{t}}+\frac{\Delta t}{k_BT}\mathbf{D\cdot F}+\sqrt{2\Delta t}\mathbf{S\cdot g}" /></a>

where **D** is the Ronte-Prager-Yamakawa diffusion tensor:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{align*}&space;\mathbf{D_{ii}}&=\frac{k_BT}{6\pi&space;\eta&space;\sigma}\mathbf{I}\\&space;\mathbf{D_{ij}}&=\frac{k_BT}{8\pi\eta}\cdot\frac{1}{r_{ij}}\cdot[(1&plus;\frac{2\sigma^2}{3r_{ij}^2})\mathbf{I}&plus;(1-\frac{2\sigma^2}{r_{ij}^2})\frac{\mathbf{r_{ij}r_{ij}}}{r_{ij}^2}]~~~~\mathrm{for}~r_{ij}\geq2\sigma\\&space;\mathbf{D_{ij}}&=\frac{k_BT}{8\pi\eta}&space;\cdot&space;\frac{1}{r_{ij}}&space;\cdot&space;[\frac{r_{ij}}{2\sigma}(\frac{8}{3}-\frac{3r_{ij}}{4\sigma})\mathbf{I}&space;&plus;&space;(\frac{r_{ij}}{4\sigma})\frac{\mathbf{r_{ij}r_{ij}}}{r_{ij}^2}]~~~~\mathrm{for}~r_{ij}<2\sigma&space;\end{align*}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\begin{align*}&space;\mathbf{D_{ii}}&=\frac{k_BT}{6\pi&space;\eta&space;\sigma}\mathbf{I}\\&space;\mathbf{D_{ij}}&=\frac{k_BT}{8\pi\eta}\cdot\frac{1}{r_{ij}}\cdot[(1&plus;\frac{2\sigma^2}{3r_{ij}^2})\mathbf{I}&plus;(1-\frac{2\sigma^2}{r_{ij}^2})\frac{\mathbf{r_{ij}r_{ij}}}{r_{ij}^2}]~~~~\mathrm{for}~r_{ij}\geq2\sigma\\&space;\mathbf{D_{ij}}&=\frac{k_BT}{8\pi\eta}&space;\cdot&space;\frac{1}{r_{ij}}&space;\cdot&space;[\frac{r_{ij}}{2\sigma}(\frac{8}{3}-\frac{3r_{ij}}{4\sigma})\mathbf{I}&space;&plus;&space;(\frac{r_{ij}}{4\sigma})\frac{\mathbf{r_{ij}r_{ij}}}{r_{ij}^2}]~~~~\mathrm{for}~r_{ij}<2\sigma&space;\end{align*}" title="\begin{align*} \mathbf{D_{ii}}&=\frac{k_BT}{6\pi \eta \sigma}\mathbf{I}\\ \mathbf{D_{ij}}&=\frac{k_BT}{8\pi\eta}\cdot\frac{1}{r_{ij}}\cdot[(1+\frac{2\sigma^2}{3r_{ij}^2})\mathbf{I}+(1-\frac{2\sigma^2}{r_{ij}^2})\frac{\mathbf{r_{ij}r_{ij}}}{r_{ij}^2}]~~~~\mathrm{for}~r_{ij}\geq2\sigma\\ \mathbf{D_{ij}}&=\frac{k_BT}{8\pi\eta} \cdot \frac{1}{r_{ij}} \cdot [\frac{r_{ij}}{2\sigma}(\frac{8}{3}-\frac{3r_{ij}}{4\sigma})\mathbf{I} + (\frac{r_{ij}}{4\sigma})\frac{\mathbf{r_{ij}r_{ij}}}{r_{ij}^2}]~~~~\mathrm{for}~r_{ij}<2\sigma \end{align*}" /></a>

**g** is a Gaussian distributed array and **S** is obtained by Cholesky decomposition of **D** since:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{D}=\left[\begin{matrix}\mathbf{D_{11}}&space;&&space;\mathbf{D_{12}}&space;&&space;\dots&space;&&space;\mathbf{D_{1n}}\\&space;\mathbf{D_{21}}&space;&&space;\mathbf{D_{22}}&space;&&space;\dots&space;&&space;\mathbf{D_{2n}}\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\mathbf{D_{n1}}&space;&&space;\mathbf{D_{n2}}&space;&&space;\dots&space;&&space;\mathbf{D_{nn}}\end{matrix}\right]=\mathbf{S^T}\mathbf{S}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\mathbf{D}=\left[\begin{matrix}\mathbf{D_{11}}&space;&&space;\mathbf{D_{12}}&space;&&space;\dots&space;&&space;\mathbf{D_{1n}}\\&space;\mathbf{D_{21}}&space;&&space;\mathbf{D_{22}}&space;&&space;\dots&space;&&space;\mathbf{D_{2n}}\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\mathbf{D_{n1}}&space;&&space;\mathbf{D_{n2}}&space;&&space;\dots&space;&&space;\mathbf{D_{nn}}\end{matrix}\right]=\mathbf{S^T}\mathbf{S}" title="\mathbf{D}=\left[\begin{matrix}\mathbf{D_{11}} & \mathbf{D_{12}} & \dots & \mathbf{D_{1n}}\\ \mathbf{D_{21}} & \mathbf{D_{22}} & \dots & \mathbf{D_{2n}}\\ \vdots & \vdots & \ddots & \vdots \\ \mathbf{D_{n1}} & \mathbf{D_{n2}} & \dots & \mathbf{D_{nn}}\end{matrix}\right]=\mathbf{S^T}\mathbf{S}" /></a>


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
