#include <iostream>
#include <iomanip>
#include <map>
#include <fstream>
#include <cmath>
#include "../include/dyna.hpp"
#include "../include/define.hpp"
#include "../include/util.hpp"

void DynaSystem::setup_dyna(std::vector<std::string> cmds, UserVar& user_var) {
    std::cout << "DynaSetup> Preparing for dynamics simulation" << std::endl;
    std::map<std::string,std::string> dyna_opt;
    for(auto it_cmd=cmds.begin()+1; it_cmd!=cmds.end(); it_cmd+=2) {
        auto retrieved_var = user_var.query(*(it_cmd+1));
        dyna_opt[*it_cmd] = retrieved_var;
    }
    std::cout << "DynaSetup> The following parameters will be used:" << std::endl;
    for(auto it_dyna=dyna_opt.begin(); it_dyna!=dyna_opt.end(); it_dyna++) {
        std::cout << "DynaSetup> " << std::left
                  << std::setw(10) << it_dyna->first  << " = " 
                  << std::setw(10) << it_dyna->second << std::endl;
    }
    /**************** For dynamics ***********************/
    this->_dyna_config.tstep   = std::stof(dyna_opt["tstep"]);
    this->_dyna_config.zeta    = std::stof(dyna_opt["zeta"]);
    this->_dyna_config.temp    = std::stof(dyna_opt["temp"]);
    this->_dyna_config.nstep   = std::stol(dyna_opt["nstep"]);
    this->_dyna_config.outfreq = std::stol(dyna_opt["outfreq"]);
    this->_dyna_config.dcdfreq = std::stol(dyna_opt["dcdfreq"]);
    this->_dyna_config.nbdfreq = std::stol(dyna_opt["nbdfreq"]);
    this->_dyna_config.dijfreq = std::stol(dyna_opt["dijfreq"]);
    this->_dyna_config.dcdname = dyna_opt["dcdname"];
    /*****************************************************/
    const size_t NATOM = this->get_NATOM();
    /**************** For dcd header *********************/
    this->_dcd_header.natom = NATOM;
    this->_dcd_header.nfile = _dyna_config.nstep/_dyna_config.dcdfreq;
    this->_dcd_header.npriv = 0;
    this->_dcd_header.nsavc = _dyna_config.dcdfreq;
    this->_dcd_header.nstep = _dyna_config.nstep;
    this->_dcd_header.ifpbc = 0;
    this->_dcd_header.delta = _dyna_config.tstep * TIMFAC;
    this->_dcd_header.prog_title = "REMARK CREATED BY EZAFM FILENAME";
    this->_dcd_header.user_title = "REMARK CREATED BY EZAFM USER";
    /*****************************************************/
    /* initialize random number array and run dynamics. */
    this->init_rand(NATOM);
    this->init_energy(NATOM);
    //if(qhydro) init_tensor(NATOM);
}

void DynaSystem::run_dyna() {
    std::cout << "DynaRun> Running dynamics for (" << _dyna_config.nstep 
              << ") steps ..." << std::endl; 
    /* coefficients for dispalcement calculation*/
    const size_t NATOM = this->get_NATOM();
    const double DELTA = _dyna_config.tstep * TIMFAC;
    const double kT = KB * _dyna_config.temp;
    const double ZETA = _dyna_config.zeta;
    const double VAR = sqrt(2*kT*ZETA/DELTA);
    const double COEFF = DELTA / ZETA;
    const double RCUT  = 16.0;
    /* initialize trajectory writing */
    std::ofstream dcd_file(_dyna_config.dcdname, std::ios::binary);
    this->write_dcdheader(dcd_file);
    /* run dynamics */
    size_t istep = 1;
    this->update_nonbond(RCUT);
    this->print_energy(istep, 1==istep);
    
    while(istep <= _dyna_config.nstep) {
        
        /* 0. update nonbonded interaction pair list */
        if(0==istep%_dyna_config.nbdfreq) {
            auto n_pairs = this->update_nonbond(RCUT);
            if(0==istep%_dyna_config.outfreq) {
                std::cout << "NonBond> Nonbonded pair list update (" 
                          << n_pairs << ") found in step (" 
                          << istep << ")" <<std::endl;
            }
                
        }
        
        /* 1. fill the random number array */
        this->gen_rand(0.0, VAR);
        
        /* 2. calculate gradients */
        this->compute_energy();
        
        /* 3. if AFM: internal energy must be calculated before this step */
        if(_afm_config.do_afm) {
            size_t index_n = _afm_config.nterm-1;
            size_t index_c = _afm_config.cterm-1;
            auto afm_coors = AfmPair {get_atom_coor(index_n), get_atom_coor(index_c)};
            auto afm_forces = apply_afm(_dyna_config.tstep, afm_coors);
            auto afm_atoms = std::vector<size_t> {index_n, index_c};
            this->other_forces<std::array<Vec3d,2> >(afm_atoms, afm_forces);
        }
        
        /* 4. At this point all force calculations are done. */
        if(0==istep%_dyna_config.outfreq)
            this->print_energy(istep, 1==istep);
        
        /* 5. calculate diffusion tensor (if necessary) */
        // if(0==istep%dijfreq)
        
        /* 6. apply displacement propagation */
        for(size_t iatom=0; iatom < NATOM; ++iatom) {
            auto tmp_rand  = get_rand(iatom);
            auto tmp_force = get_force(iatom);
            auto dr = COEFF*(tmp_force + tmp_rand);
            if(is_movable(iatom))
                this->move_cor(dr, iatom);
        }
        
        /* 7. write coordinate to trajectory */
        if(0==istep%_dyna_config.dcdfreq)
            this->write_dcdframe(dcd_file, p_xcoor(), p_ycoor(), p_zcoor(), NATOM);
        istep++;
    }
    dcd_file.close();
}