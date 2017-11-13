#ifndef USERVAR_H
#define USERVAR_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "../include/main.h"

//default req'd variables
const struct Default_Var {
    real tstep = 1.0;
    real zeta = 50.0;
    real temp = 298.0;
    unsigned int nstep = 1000;
    unsigned int iseed = 1234;
    unsigned int outfreq = 10;
    unsigned int dcdfreq = 10;
    unsigned int nbdfreq = 10;
    unsigned int dijfreq = 10;
    std::string hydro = "none";
} default_var;

void set_uservar(std::vector<std::string> cmds);

std::string get_uservar(std::string name);

void print_uservar();

#endif
