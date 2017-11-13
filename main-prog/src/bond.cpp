#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../include/main.h"
#include "../include/energy.h"
#include "../include/bond.h"
#include "../include/param.h"
#include "../include/psf.h"
#include "../include/coor.h"

real bond_energy () {
    real e_bond = 0.0;
    if(psf_bond().empty()) return e_bond;
    for(auto it=psf_bond().cbegin(); it!=psf_bond().cend(); ++it) {
        /* retrieve energy parameters */
        unsigned int atom_i = it->atom_i;
        unsigned int atom_j = it->atom_j;        
        std::string  type_i = get_atom_type(atom_i);
        std::string  type_j = get_atom_type(atom_j);
        real k_bond = get_bond_params(PrmBondType(type_i,type_j)).k_bond;
        real r_eq   = get_bond_params(PrmBondType(type_i,type_j)).r_eq;
        /* retrieven bond component coordinates */
        auto i = atom_i-1; auto j = atom_j-1;
        Vector r_i(get_xcoor(i),get_ycoor(i),get_zcoor(i));
        Vector r_j(get_xcoor(j),get_ycoor(j),get_zcoor(j));
        Vector r_ij = r_i - r_j;
        real r_diff = r_ij.norm() - r_eq;
        /* calculate force and energy */
        real force  = k_bond * r_diff;
        real energy = force  * r_diff;
        Vector f_i = -2.0f * force * r_ij.unitvec(); 
        Vector f_j = -1.0f * f_i;
        /* add forces to gradient array */
        add_to_force(i, f_i);
        add_to_force(j, f_j);
        e_bond += energy;
    }
    return e_bond;
}