#ifndef ENERGY_H
#define ENERGY_H

#include "../include/main.h"

struct Energy_Info {
    real e_bond;
    real e_angl;
    real e_dihe;
    real e_nbnd;
    real e_total;
};

const Energy_Info& get_energy_info();

bool is_ready_for_energy();

void init_energy(unsigned int natom);

void calc_energy();

const real& get_force(unsigned int query);

real  get_norm_grad();

void add_to_force(const unsigned int& index, const Vector& force);

#endif