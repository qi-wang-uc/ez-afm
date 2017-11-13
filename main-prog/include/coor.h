#ifndef COOR_H
#define COOR_H

#include <string>
#include "../include/main.h"

bool read_cor(const std::string &inp_name);

// In this case, a good COR file should have at least one atom
// TODO: what is a good COR file?
bool is_good_cor();

unsigned int get_natom_cor();

void print_coor();

real* p_xcoor();
real* p_ycoor();
real* p_zcoor();

const real& get_xcoor(const unsigned int atomid);
const real& get_ycoor(const unsigned int atomid);
const real& get_zcoor(const unsigned int atomid);

void add_to_coor(const unsigned int atomid, const Vector& coor);

#endif
