#ifndef DYNA_H
#define DYNA_H

#include <vector>
#include <string>

/*
struct Dyna_Args {
    real tstep; // timestep
    real zeta;
    real temp;            
    unsigned int nstep; // should be converted from <real>. 
    unsigned int outfreq; 
    unsigned int dcdfreq; 
    unsigned int nbdfreq; 
    unsigned int dijfreq;
    unsigned int resisep;
    int qhydro; // enum?
    std::string dcdname;
};
*/

bool prep_dyna(std::vector<std::string> cmds);

void run_dyna(unsigned int nstep, real tstep, real zeta, real temp, 
               unsigned int outfreq, unsigned int dcdfreq, unsigned int nbdfreq, 
               unsigned int dijfreq, int qhydro, std::string dcdname);

#endif            