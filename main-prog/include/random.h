#ifndef RAND_H
#define RAND_H

#include "../include/main.h"

void gen_random(real mean, real dev);

void init_random(const unsigned int natom);

const real& get_random(unsigned int query);

void print_random();

#endif