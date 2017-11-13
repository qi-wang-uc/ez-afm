#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../include/main.h"
#include "../include/afm.h"
#include "../include/psf.h"
#include "../include/dcd.h"
#include "../include/util.h"
#include "../include/coor.h"
#include "../include/dyna.h"
#include "../include/nonbond.h"
#include "../include/energy.h"
#include "../include/random.h"
#include "../include/uservar.h"

bool prep_dyna(std::vector<std::string> cmds) {
    std::cout << "PREPDYNA> Preparing for dynamics simulation." << std::endl;
    std::map<std::string,std::string> dyna_opt;
    for(auto it_cmd=cmds.begin()+1; it_cmd!=cmds.end(); it_cmd+=2) {
        auto retrieved_var = get_uservar(*(it_cmd+1));
        if(retrieved_var.empty()) return false;
        dyna_opt[*it_cmd] = retrieved_var;
    }
    std::cout << "PREPDYNA> The following parameters will be used for dynamics simulations:" << std::endl;
    std::cout << "******************************" << std::endl;
    for(auto it_dyna=dyna_opt.begin(); it_dyna!=dyna_opt.end(); it_dyna++) {
            std::cout << std::left << std::setw(10) <<"[" + it_dyna->first + "]" << " is set to " 
                 << std::setw(10) << it_dyna->second << std::endl;
    }
    std::cout << "******************************" << std::endl;
    //TODO: cout << "xxx not found in input file, yyy will be used."
    // Now convert the strings to functional variable types.
    int qhydro = 0; // default: when no hydrodynamics interactions.
    if(dyna_opt["hydro"]=="cholesky") qhydro = 1;
    if(dyna_opt["hydro"]=="tea") qhydro = 2;
    unsigned int nstep = std::stol(dyna_opt["nstep"]);
    //unsigned int iseed = std::stol(dyna_opt["iseed"]);
    unsigned int outfreq = std::stol(dyna_opt["outfreq"]);
    unsigned int dcdfreq = std::stol(dyna_opt["dcdfreq"]);
    unsigned int nbdfreq = std::stol(dyna_opt["nbdfreq"]);
    unsigned int dijfreq = std::stol(dyna_opt["dijfreq"]);
    real tstep = std::stof(dyna_opt["tstep"]);
    real zeta = std::stof(dyna_opt["zeta"]);
    real temp = std::stof(dyna_opt["temp"]);
    std::string dcdname = dyna_opt["dcdname"];
    /* initialize random number array and run dynamics. */
    init_random(get_psf_info().natom);
    init_energy(get_psf_info().natom);
    const unsigned int resi_sep = 2;  // residue index seperation.
    init_pairlist(get_psf_info().natom, resi_sep);
    //print_grad();
    //print_coor();
    //if(qhydro) init_tensor(get_psf_info().natom);
    run_dyna(nstep, tstep, zeta, temp, outfreq, dcdfreq, nbdfreq, dijfreq, qhydro, dcdname);
    return true;
}

void run_dyna(unsigned int nstep, real tstep, real zeta, real temp, unsigned int outfreq, 
    unsigned int dcdfreq, unsigned int nbdfreq, unsigned int dijfreq, int qhydro,std::string dcdname) {
    std::cout << "RUNDYNA> Running dynamics for (" << nstep << ") steps ..." << std::endl; 
    const unsigned int resi_sep = 2;  // residue index seperation.
    /* coefficients for dispalcement calculation*/
    const real kB = 0.00198993;
    const real timfac = 4.88882129e-02; // convert AKMA time unit to picosecond.
    const real delta = tstep * timfac;
    const real kT = kB * temp;
    const real coeff = delta / zeta;
    const real cutoff = 16.0;   // cutoff value for building non-bonded pair list.
    /* initialize trajectory writing */
    DCD_Header dcd_header = {
        get_psf_info().natom,   //NATOM
        nstep/dcdfreq,          //NFILE
        0,                      //NPRIV
        dcdfreq,                //NSAVC
        nstep,                  //NSTEP
        0,                      //IF PBC
        delta,                   //DELTA
        "REMARK CREATED BY NAMD FILENAME",       //PROGRAM TITLE
        "REMARK CREATED BY NAMD USER",           //USER TITLE
    };
    DCD_Pads dcd_pads;
    std::ofstream dcd_file(dcdname, std::ios::binary);
    write_dcdheader(dcd_file, dcd_header, dcd_pads);
    /* run dynamics */
    unsigned int istep = 1;
    unsigned int natom = get_psf_info().natom;
    const real var = sqrt(2*kT*zeta/delta);
    update_pairlist(cutoff, resi_sep);
    std::cout << "ENEINFO> " << std::setw(10) << "STEP" << std::setw(16) << "EBOND" 
        << std::setw(16) << "E_ANGLE" << std::setw(16) << "E_DIHEDRAL" 
        << std::setw(16) << "E_NBOND" << std::setw(16) <<"E_TOTAL" << std::endl;
    /* Dynamics start from here */
    while(istep <= nstep) {
        /* 0. update nonbonded interaction pair list */
        if(0==istep%nbdfreq) {
            auto n_pairs = update_pairlist(cutoff, resi_sep);
            std::cout << "NONBOND> Nonbonded pair list update (" << n_pairs << ") found in step (" 
                      << istep << ")" <<std::endl;
        }
        /* 1. fill the random number array */
        gen_random(0.0, var);
        /* 2. calculate gradients */
        calc_energy();
        /* 3. if do AFM. Note: internal energy must be calculated first before this step because
            energy array is zeroed out in calc_energy(). */
        if(get_afm_info().do_afm)
            apply_afm(tstep);

        if(0==istep%outfreq) { 
            std::cout << "ENEINFO> " 
                << std::setw(10) << istep 
                << std::setw(16) << get_energy_info().e_bond 
                << std::setw(16) << get_energy_info().e_angl 
                << std::setw(16) << get_energy_info().e_dihe 
                << std::setw(16) << get_energy_info().e_nbnd 
                << std::setw(16) << get_energy_info().e_total 
                << std::endl;
        }
        
        /* 4. calculate diffusion tensor (if necessary) */
        // if(0==istep%dijfreq)

        /* 5. apply displacement propagation */
        for(unsigned int iatom=0; iatom < natom; ++iatom) {
            Vector incre(coeff * (get_force(3*iatom+0) + get_random(3*iatom+0)),
                         coeff * (get_force(3*iatom+1) + get_random(3*iatom+1)),
                         coeff * (get_force(3*iatom+2) + get_random(3*iatom+2)));
            add_to_coor(iatom, incre);
        }

        /* 6. write coordinate to trajectory */
        if(0==istep%dcdfreq)
            write_dcdframe(dcd_file, natom, dcd_pads);
        
        istep++;
    }
    dcd_file.close();
}
