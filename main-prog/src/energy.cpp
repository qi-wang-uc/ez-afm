#include <iostream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include "../include/param.h"
#include "../include/psf.h"
#include "../include/coor.h"
#include "../include/energy.h"
#include "../include/bond.h"
#include "../include/angle.h"
#include "../include/dihedral.h"
#include "../include/nonbond.h"

/* Allocate memory space for gradients (forces).
   Then calculate each energy terms. */
static std::vector<real> grad;    //size = 3N
static Energy_Info energy_info;

void init_energy(unsigned int natom) {
    if((!is_ready_for_energy())||(!grad.empty())) return;
    grad.resize(3*natom);
}

void calc_energy() {
    /* zero out force array */
    std::fill(grad.begin(),grad.end(), 0.0);
    /* calculate each energy term */
    energy_info.e_bond = bond_energy();
    //energy_info.e_angl = angl_energy();
    //energy_info.e_dihe = dihe_energy();
    energy_info.e_nbnd = nbnd_energy();
    energy_info.e_total = energy_info.e_bond + 
                          energy_info.e_angl + 
                          energy_info.e_dihe +
                          energy_info.e_nbnd;
}

const real& get_force(unsigned int query) {
    return grad[query];
}

real get_norm_grad() {
    real dp = std::inner_product(grad.cbegin(),grad.cend(),grad.cbegin(),0.0);
    return dp / grad.size();
}

bool is_ready_for_energy() {
    if(!is_good_params()) return false;
    if(!is_good_psf()) return false;
    if(!is_good_cor()) return false;
    if(get_natom_cor()!=get_psf_info().natom) {
        std::cout << "ERROR> Atom numbers in PSF and COR files do not match." << std::endl;
        return false;
    }
    return true;
}

const Energy_Info& get_energy_info() {
    return energy_info;
}

void add_to_force(const unsigned int& atomid, const Vector& force) {
    grad[3*atomid+0] += force.x;
    grad[3*atomid+1] += force.y;
    grad[3*atomid+2] += force.z;
}