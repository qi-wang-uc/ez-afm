#ifndef NBOND_H
#define NBOND_H

#include "../include/main.h"

void init_pairlist(unsigned int natom, unsigned int resi_sep);

unsigned int update_pairlist(real cutoff, unsigned int resi_sep);

inline unsigned int get_pairlist_index(unsigned int i, unsigned int j, unsigned int s, unsigned int n);

real nbnd_energy();

#endif
