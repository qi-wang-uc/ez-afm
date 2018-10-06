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
    Int  nterm = 1;
    Int  cterm = 1;
    Real k_afm = 0;
    Real v_afm = 0;
    Real tmp_dist = 0.0;
    Real max_dist = 0.0;
};

class AFM {
    private:
        AfmConfig _afm_config;
    public:
        void setup_afm(StrVec cmds, UserVar& user_var);
        AfmPair apply_afm(const Real& tstep, const AfmPair& afm_coors);
        AfmConfig const& get_config() const;
};

#endif