#include <iostream>
#include <iomanip>
#include "../include/sd.h"
#include "../include/energy.h"
#include "../include/coor.h"

void steepest_descent(const unsigned int nstep) {
    /* Banner of energy minimization output */
    const int w = 10;
    std::cout << "EMINI>" << std::setw(w) << "STEP" << std::setw(w) 
        << "E_TOTAL" << std::setw(w) << "E_BOND" << std::setw(w) 
        << "E_ANGL" << std::setw(w) << "E_DIHE" << std::setw(w) 
        << "E_NBND" << std::endl;
    /* Implement the algorithm */
    unsigned int istep = 1;
    unsigned int natom = get_natom_cor();
    const real stepsize = 0.02;
    while(istep <= nstep) {
        /* 1. claculate the gradient. */
        calc_energy();
        /* 2. caculate normalized dot product of gradient array. */
        real norm_grad = get_norm_grad();
        // 3. apply steepest descent algorithm.
        for(unsigned int i=0; i<natom; ++i) {
            /*
            get_xcoor()[i] -= stepsize * get_grad()[i+0] / norm_grad;
            get_ycoor()[i] -= stepsize * get_grad()[i+1] / norm_grad;
            get_zcoor()[i] -= stepsize * get_grad()[i+2] / norm_grad;
            */
        }
        // 4. print out energy of this step.
        std::cout << "EMINI>" << std::setw(w) << istep << std::setw(w) 
            << get_energy_info().e_total << std::setw(w) 
            << get_energy_info().e_bond << std::setw(w) 
            << get_energy_info().e_angl << std::setw(w) 
            << get_energy_info().e_dihe << std::setw(w) 
            << get_energy_info().e_nbnd << std::endl;
        istep++;
    }
}