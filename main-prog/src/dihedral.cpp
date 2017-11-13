#include <iostream>
#include <string>
#include <cmath>
#include "../include/main.h"
#include "../include/energy.h"
#include "../include/dihedral.h"
#include "../include/param.h"
#include "../include/psf.h"
#include "../include/coor.h"
#include "../include/util.h"

/* TODO: not numerically stable, need to improve algorithm
*/

real dihe_energy () {
    real e_dihe = 0.0;
    if(psf_dihe().empty()) return e_dihe;
    for(auto it=psf_dihe().cbegin(); it!=psf_dihe().cend(); ++it) {
        /* retrieve atom information */
        unsigned int atom_i = it->atom_i;
        unsigned int atom_j = it->atom_j;
        unsigned int atom_k = it->atom_k;
        unsigned int atom_l = it->atom_l;
        std::string  type_i = get_atom_type(atom_i);
        std::string  type_j = get_atom_type(atom_j);
        std::string  type_k = get_atom_type(atom_k);
        std::string  type_l = get_atom_type(atom_l);
        auto i=atom_i-1; auto j=atom_j-1; auto k=atom_k-1; auto l=atom_l-1;
        // TODO: use const version of coordinate getter.
        Vector r_i(get_xcoor(i), get_ycoor(i), get_zcoor(i));
        Vector r_j(get_xcoor(j), get_ycoor(j), get_zcoor(j));
        Vector r_k(get_xcoor(k), get_ycoor(k), get_zcoor(k));
        Vector r_l(get_xcoor(l), get_ycoor(l), get_zcoor(l));
        auto r_ij = r_i - r_j; 
        auto r_jk = r_j - r_k; auto r_kj = r_k - r_j; 
        auto r_kl = r_k - r_l; auto r_lk = r_l - r_k;
        real rjk2inv = 1.0f / r_jk.norm2();
        /* vectors R and S form the "eigen-angle" of dihedral such that cos(<R,S>)=cos(phi)*/
        Vector R = r_ij - (r_ij.dot_product(r_kj.unitvec()))*r_kj.unitvec();
        Vector S = r_lk - (r_lk.dot_product(r_kj.unitvec()))*r_kj.unitvec();
        real cos_phi = R.unitvec().dot_product(S.unitvec());
        /* vectors M and N are used to determine the sign of phi using IUPAC convention */
        Vector M = r_ij.cross_product(r_kj);
        Vector N = r_kj.cross_product(r_kl);
        real dp = r_kj.dot_product(M.cross_product(N));
        real sign_of_phi = dp > 0.0 ? 1.0 : -1.0;
        real phi = sign_of_phi * acos(cos_phi);
        /* Since f_i+f_j+f_k+f_l=0, choose vector T such that f_j=-f_i+T and f_k=-f_l-T. */
        auto dihe_prm = get_dihe_params(PrmDiheType(type_i,type_j,type_k,type_l));
        for(auto itd=dihe_prm.cbegin(); itd!=dihe_prm.cend(); ++itd) {
            real k_dihe = itd->k_dihe;
            unsigned short mul = itd->mul;
            real energy = k_dihe * (1.0f + cos(mul*phi));
            real force = -1.0f * k_dihe * mul * sin(mul*phi);
            Vector f_i = force * 1.0f/R.norm() * (S.unitvec()-cos_phi*R.unitvec());
            Vector f_l = force * 1.0f/S.norm() * (R.unitvec()-cos_phi*S.unitvec());
            Vector T = rjk2inv*(r_ij.dot_product(r_kj)*f_i-(r_kl.dot_product(r_kj)*f_l));
            Vector f_j = -1.0f * f_i + T;
            Vector f_k = -1.0f * f_l - T;
            /* apply changes to force array */
            add_to_force(i, f_i);
            add_to_force(j, f_j);
            add_to_force(k, f_k);
            add_to_force(l, f_l);
            e_dihe += energy;
        }
    }
    return e_dihe;
}











