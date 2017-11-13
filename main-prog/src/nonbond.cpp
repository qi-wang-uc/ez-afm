#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <algorithm>
#include "../include/nonbond.h"
#include "../include/param.h"
#include "../include/psf.h"
#include "../include/coor.h"
#include "../include/energy.h"

//static std::map<std::pair<unsigned int, unsigned int>, unsigned short> pair_list;


static std::vector<unsigned short> pair_list;

void init_pairlist(unsigned int natom, unsigned int resi_sep) {
    auto n_possible_pairs = (natom-resi_sep)*(natom-resi_sep+1)/2;
    pair_list.resize(n_possible_pairs);
    std::cout << "After initialization, (" << n_possible_pairs 
              << ") possible non-bonded pairs were found." << std::endl;
}

unsigned int update_pairlist(real cutoff, unsigned int resi_sep) {
    std::fill(pair_list.begin(), pair_list.end(), 0);
    real cutoff2 = cutoff*cutoff;
    unsigned int pair_list_counter = 0;
    const unsigned int natom = get_psf_info().natom;
    const unsigned int natom_sep = natom - resi_sep;
    for(unsigned int i=0; i < natom-resi_sep; i++) {
        Vector r_i(get_xcoor(i),get_ycoor(i),get_zcoor(i));
        for(unsigned int j=i+resi_sep; j < natom; j++) {
       //     std::cout << "DEBUG >>>> " << i << " - " << j << " - " << resi_sep << " - " << natom << " : " << get_pairlist_index(i,j,resi_sep,natom_sep) << std::endl;
            Vector r_j(get_xcoor(j),get_ycoor(j),get_zcoor(j));
            Vector r_ij = r_i - r_j;
            if (r_ij.norm2() < cutoff2) {
                //pair_list.at(get_pairlist_index(i,j,resi_sep,natom_sep)) = 1;
                pair_list[get_pairlist_index(i,j,resi_sep,natom_sep)] = 1;
                pair_list_counter++;
            }
        }
    }
    return pair_list_counter;
}

real nbnd_energy () {
    real e_nbnd = 0.0;
    if(pair_list.empty()) return e_nbnd;
    const unsigned int natom = get_psf_info().natom;
    const unsigned int resi_sep = 2;
    const unsigned int natom_sep = natom - resi_sep;
    const real sigma   = 3.8;
    const real epsilon = 2.0;
    real sigma2 = sigma * sigma;

    for(unsigned int i=0; i<natom; ++i) {
        for(unsigned int j=i+resi_sep; j<natom-resi_sep; ++j) {
            if(1==get_pairlist_index(i,j,resi_sep,natom_sep)) {
                Vector r_i(get_xcoor(i), get_ycoor(i), get_zcoor(i));
                Vector r_j(get_xcoor(j), get_ycoor(j), get_zcoor(j));
                Vector r_ij = r_i - r_j;
                real fr2 = sigma2 / r_ij.norm2();
                real fr6 = fr2 * fr2 * fr2;
                real energy = 4.0*epsilon*fr6*(fr6-1.0);
                real force = 48.0*epsilon*fr6*(fr6-0.5)/r_ij.norm2();
                Vector f_i = force * r_ij;
                Vector f_j = -1.0f * f_i;
                // apply changes into force array
                add_to_force(i, f_i);
                add_to_force(j, f_j);
                e_nbnd += energy;
            }
        }
    } 
    return e_nbnd;
}

unsigned int get_pairlist_index(unsigned int i, unsigned int j, unsigned int s, unsigned int n) {
	return (2*n-i+1)*i/2 + (j-(i+s)+1) - 1;
}

/*
*/


    /*
    for(auto it=pair_list.cbegin(); it!=pair_list.cend(); ++it) {
        if(0==it->second) continue;
        auto atom_i = it->first.first;
        auto atom_j = it->first.second;
        auto i = atom_i - 1;
        auto j = atom_j - 1;
        Vector r_i(get_xcoor(i), get_ycoor(i), get_zcoor(i));
        Vector r_j(get_xcoor(j), get_ycoor(j), get_zcoor(j));
        Vector r_ij = r_i - r_j;
        real fr2 = sigma2 / r_ij.norm2();
        real fr6 = fr2 * fr2 * fr2;
        real energy = 4.0*epsilon*fr6*(fr6-1.0);
        real force = 48.0*epsilon*fr6*(fr6-0.5)/r_ij.norm2();
        Vector f_i = force * r_ij;
        Vector f_j = -1.0f * f_i;
    
        // apply changes into force array
        add_to_force(i, f_i);
        add_to_force(j, f_j);
        e_nbnd += energy;
    }*/