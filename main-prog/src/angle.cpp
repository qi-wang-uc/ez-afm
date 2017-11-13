#include <iostream>
#include <iomanip>
#include <cmath>
#include "../include/main.h"
#include "../include/energy.h"
#include "../include/angle.h"
#include "../include/param.h"
#include "../include/psf.h"
#include "../include/coor.h"

/* TODO: not numerically stable, need to improve algorithm */

real angl_energy () {
    real e_angle = 0.0;
    if(psf_angl().empty()) return e_angle;
    for(auto it=psf_angl().cbegin(); it!=psf_angl().cend(); ++it) {
        /* retrieve angle energy parameters */
        unsigned int atom_i = it->atom_i;
        unsigned int atom_j = it->atom_j;
        unsigned int atom_k = it->atom_k;
        std::string  type_i = get_atom_type(atom_i);
        std::string  type_j = get_atom_type(atom_j);
        std::string  type_k = get_atom_type(atom_k);
        real k_angle  = get_angl_params(PrmAnglType(type_i,type_j,type_k)).k_angle;
        real theta_eq = get_angl_params(PrmAnglType(type_i,type_j,type_k)).theta_eq;
        /* retrieven angle component coordinates */
        auto i = atom_i-1; auto j = atom_j-1; auto k = atom_k-1;
        Vector r_i(get_xcoor(i), get_ycoor(i), get_zcoor(i));
        Vector r_j(get_xcoor(j), get_ycoor(j), get_zcoor(j));
        Vector r_k(get_xcoor(k), get_ycoor(k), get_zcoor(k));
        Vector r_ij = r_i - r_j;
        Vector r_kj = r_k - r_j;
        real cos_theta = r_ij.dot_product(r_kj)/(r_ij.norm()*r_kj.norm());
        real sin_theta = sqrt(1-cos_theta*cos_theta);
        real theta = acos(cos_theta);
        /* calculate force and energy */
        real force  = k_angle * (theta - theta_eq);
        real energy = -1.0f * force * (theta - theta_eq);
        force /= sin_theta;
        Vector f_i = (force/r_ij.norm())*(r_kj.unitvec()-cos_theta*r_ij.unitvec());
        Vector f_k = (force/r_kj.norm())*(cos_theta*r_kj.unitvec()-r_ij.unitvec());
        Vector f_j = -1.0f * (f_i + f_k);
        /* apply changes into force array */
        add_to_force(i, f_i);
        add_to_force(j, f_j);
        add_to_force(k, f_k);
        e_angle += energy;
    }
    return e_angle;
}