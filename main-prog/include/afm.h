#ifndef AFM_H
#define AFM_H

#include <vector>
#include <string>
#include "main.h"

struct Afm_Info {
    bool do_afm = false;
    unsigned int nterm = 1;
    unsigned int cterm = 1;
    real k_afm = 100.0;
    real v_afm = 1e-3;
    real tmp_dist = 0.0;
    real max_dist = 0.0;
};

const Afm_Info& get_afm_info();

bool prep_afm(std::vector<std::string> cmds);

void apply_afm(real tstep);



#endif