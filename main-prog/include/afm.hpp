#ifndef AFM_HPP
#define AFM_HPP

#include <vector>
#include <string>
#include <array>
#include "define.hpp"
#include "uservar.hpp"
/* All external forces go here, AFM, restraints, etc.*/

// Pair-wise coodinates and forces. Do NOT use std::pair here
using AfmPair = std::array<Vec3d, 2>;

struct AfmConfig {
    bool do_afm = false;
    size_t nterm = 1;
    size_t cterm = 1;
    double k_afm = 0;
    double v_afm = 0;
    double tmp_dist = 0.0;
    double max_dist = 0.0;
};

class AFM {
    protected:
        AfmConfig _afm_config;
    public:
        void setup_afm(std::vector<std::string> cmds, UserVar& user_var);
        AfmPair apply_afm(const double& tstep, const AfmPair& afm_coors);
};

#endif