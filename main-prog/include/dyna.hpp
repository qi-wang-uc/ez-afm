#ifndef DYNA_HPP
#define DYNA_HPP

#include <vector>
#include <string>
#include "energy.hpp"
#include "rand.hpp"
#include "afm.hpp"
#include "dcd.hpp"
#include "uservar.hpp"

struct DynaConfig {
    double tstep   = 1;
    double zeta    = 50;
    double temp    = 300;            
    size_t nstep   = 1000;
    size_t outfreq = 100; 
    size_t dcdfreq = 100; 
    size_t nbdfreq = 100; 
    size_t dijfreq = 100; 
    size_t resisep = 100; 
    size_t qhydro  = 0;
    std::string dcdname = "default";
};

class DynaSystem : public Energy, public Rand, public DcdData, public AFM {
    protected:
        DynaConfig _dyna_config;
    public:
        void setup_dyna(std::vector<std::string> cmds, UserVar& user_var);
        void run_dyna();
};

#endif
