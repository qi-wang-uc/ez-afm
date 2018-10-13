#include <iostream>
#include <iomanip>
#include <map>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "../include/dyna.hpp"
#include "../include/define.hpp"
#include "../include/util.hpp"
#include "../include/dcd.hpp"

void DynaSystem::set_dynaconfig(DynaConfig&& inp_config) {
    this->_dyna_config = std::move(inp_config);
}

void DynaSystem::setup_dyna(StrVec cmds, UserVar& user_var) {
    std::cout << "DynaSetup> Preparing for dynamics simulation" << std::endl;
    std::map<std::string,std::string> dyna_opt;
    for(auto it_cmd=cmds.begin()+1; it_cmd!=cmds.end(); it_cmd+=2) {
        auto retrieved_var = user_var.query(*(it_cmd+1));
        dyna_opt[*it_cmd] = retrieved_var;
    }
    std::cout << "DynaSetup> The following parameters will be used:" << std::endl;
    for(auto it_dyna=dyna_opt.begin(); it_dyna!=dyna_opt.end(); it_dyna++) {
        std::cout << "DynaSetup> " 
                  << std::left
                  << std::setw(10) << it_dyna->first
                  << " = " 
                  << std::setw(10) << it_dyna->second
                  << std::endl;
    }
    /**************** For dynamics ***********************/
    this->set_dynaconfig(DynaConfig(
        static_cast<Real>(std::stof(dyna_opt["tstep"])),
        static_cast<Real>(std::stof(dyna_opt["zeta"])),
        static_cast<Real>(std::stof(dyna_opt["temp"])),
        static_cast<Int> (std::stol(dyna_opt["nstep"])),
        static_cast<Int> (std::stol(dyna_opt["outfreq"])),
        static_cast<Int> (std::stol(dyna_opt["dcdfreq"])),
        static_cast<Int> (std::stol(dyna_opt["nbdfreq"])),
        static_cast<Int> (std::stol(dyna_opt["dijfreq"])),
        static_cast<Int> (std::stol(dyna_opt["hydro"])),
        static_cast<Str> (dyna_opt["dcdname"])));
    /*****************************************************/
    
    /**************** For dcd header *********************/
    this->dcd.set_dcdheader(DCD_Header(
        static_cast<int32_t>(this->psf.get_NATOM()), 
        static_cast<int32_t>(this->_dyna_config.nstep/this->_dyna_config.dcdfreq), 
        static_cast<int32_t>(0), 
        static_cast<int32_t>(this->_dyna_config.dcdfreq), 
        static_cast<int32_t>(this->_dyna_config.nstep), 
        static_cast<int32_t>(0), 
        static_cast<float>(this->_dyna_config.tstep * TIMFAC),
        static_cast<const char*>(("REMARK CREATED AT " + time_stamp()).c_str()), 
        static_cast<const char*>(("REMARK CREATED BY " + Str(getenv("USER"))).c_str())));
    /*****************************************************/
    /* initialize random number array and run dynamics. */
    const Int NATOM = this->psf.get_NATOM();
    this->rand.init_rand(NATOM);
    this->ener.init_energy(NATOM);
    if(this->_dyna_config.hydro) this->HI.init(NATOM);
}

void DynaSystem::run_dyna() {
    std::cout << "DynaRun> Running dynamics for " 
              << encap(this->_dyna_config.nstep)
              << " steps ..."
              << std::endl; 
    /* coefficients for dispalcement calculation*/
    const Int  NATOM = this->psf.get_NATOM();
    const Real DELTA = this->_dyna_config.tstep * TIMFAC;
    const Real kT    = this->_dyna_config.temp * KB;
    const Real ZETA  = this->_dyna_config.zeta;
    const Real VAR   = sqrt(2*kT*ZETA/DELTA);
    const Real COEFF = DELTA / ZETA;
    const Real RCUT  = 16.0;
    /* initialize trajectory writing */
    std::ofstream dcd_file(this->_dyna_config.dcdname, std::ios::binary);
    this->dcd.write_dcdheader(dcd_file);
    /* run dynamics */
    Int istep = 1;
    this->ener.update_nonbond(RCUT, this->cor);
    this->ener.print_energy(istep, 1==istep);
    while(istep <= this->_dyna_config.nstep) {
        /* 0. update nonbonded interaction pair list */
        if(0==istep%this->_dyna_config.nbdfreq) {
            auto n_pairs = this->ener.update_nonbond(RCUT, this->cor);
            print_nonbond(istep, n_pairs, this->_dyna_config.outfreq);    
        }
        /* 1. fill the random number array */
        this->rand.gen_rand(0.0, VAR);
        /* 2. calculate gradients */
        this->ener.compute_energy(this->psf, this->prm, this->cor);
        /* 3. if AFM: internal energy must be calculated before this step */
        if(this->afm.get_config().do_afm) {
            Int index_n = this->afm.get_config().nterm-1;
            Int index_c = this->afm.get_config().cterm-1;
            auto afm_coors = AfmPair {this->cor.get_atom_coor(index_n), this->cor.get_atom_coor(index_c)};
            auto afm_forces = this->afm.apply_afm(_dyna_config.tstep, afm_coors);
            auto afm_atoms = std::vector<Int> {index_n, index_c};
            this->ener.other_forces<std::array<Vec3d,2> >(afm_atoms, afm_forces);
        }
        /* 4. At this point all force calculations are done. */
        if(0==istep%this->_dyna_config.outfreq)
            this->ener.print_energy(istep, 1==istep);
        /* 5. Calculate displacement */
        if(1==this->_dyna_config.hydro) {
            if((0==istep%this->_dyna_config.dijfreq) || 1==istep) {
                HI.build(this->cor, this->_dyna_config.zeta);
                HI.cholesky();
            }
            HI.apply_disp_d(DELTA, this->psf, this->ener, this->cor);
            HI.apply_disp_r(DELTA, this->psf, this->ener, this->cor, this->rand);
        } else if (0==this->_dyna_config.hydro) {
            for(Int iatom=0; iatom < NATOM; ++iatom) {
                auto tmp_rand  = this->rand.get_rand(iatom);
                auto tmp_force = this->ener.get_force(iatom);
                auto dr = COEFF*(tmp_force + tmp_rand);
                this->cor.move_cor(dr, iatom, psf.is_movable(iatom));
            }
        } else {
            std::cout << "ERROR> Unknown dynamics configuration: " 
                      << encap(this->_dyna_config.hydro)
                      << std::endl;
            break;
        }
        /* 7. write coordinate to trajectory */
        if(0==istep%_dyna_config.dcdfreq || 1==istep)
            this->dcd.write_dcdframe(dcd_file, this->cor.px(), this->cor.py(), this->cor.pz(), NATOM);
        istep++;
    }
    dcd_file.close();
}

void print_nonbond(Int istep, Int n_pairs, Int out_freq) {
    if(0==istep%out_freq) {
        std::cout << "NonBond> Nonbonded pair list update " 
                  << encap(n_pairs)
                  << " found in step " 
                  << encap(istep)
                  << std::endl;
    }
}