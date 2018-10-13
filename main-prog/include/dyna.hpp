#ifndef DYNA_HPP
#define DYNA_HPP

#include <vector>
#include <string>
#include "define.hpp"
#include "energy.hpp"
#include "rand.hpp"
#include "afm.hpp"
#include "dcd.hpp"
#include "uservar.hpp"
#include "tensor.hpp"

struct DynaConfig {
    Real tstep   = 1;
    Real zeta    = 50;
    Real temp    = 300;            
    Int  nstep   = 1000;
    Int  outfreq = 100; 
    Int  dcdfreq = 100;
    Int  nbdfreq = 100; 
    Int  dijfreq = 100; 
    Int  hydro  = 0;
    Str  dcdname = "";
    DynaConfig() {}
    DynaConfig(Real tstep, Real zeta, Real temp, 
            Int nstep, Int outfreq, Int dcdfreq, Int nbdfreq, Int dijfreq, Int hydro,
            Str dcdname) :
        tstep(tstep), zeta(zeta), temp(temp),nstep(nstep),
            outfreq(outfreq), dcdfreq(dcdfreq), nbdfreq(nbdfreq), dijfreq(dijfreq), hydro(hydro),
            dcdname(dcdname) {}
};

class DynaSystem {
    private:
        DynaConfig _dyna_config;
    public:
        PsfData psf;
        PrmData prm;
        CorData cor;
        DcdData dcd;
        Energy ener;
        Rand   rand;
        AFM     afm;
        HItensor HI;
        void set_dynaconfig(DynaConfig&& inp_config);
        void setup_dyna(StrVec cmds, UserVar& user_var);
        void run_dyna();
};

void print_nonbond(Int istep, Int nonbonds, Int out_freq);

#endif
