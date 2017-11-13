#include <iostream>
#include <iomanip>
#include <map>
#include "../include/energy.h"
#include "../include/main.h"
#include "../include/coor.h"
#include "../include/afm.h"
#include "../include/psf.h"
#include "../include/uservar.h"

static Afm_Info afm_info;
/* Converts forces from (kcal/mol)/A to pN  */
const real fconv = 69.478508;

bool prep_afm(std::vector<std::string> cmds) {
    std::cout << "PREPAFM> Preparing for AFM-like simulation setup." << std::endl;
    std::map<std::string,std::string> afm_opt;
    for(auto it_cmd=cmds.begin()+1; it_cmd!=cmds.end(); it_cmd+=2) {
        auto retrieved_var = get_uservar(*(it_cmd+1));
        if(retrieved_var.empty()) return false;
        afm_opt[*it_cmd] = retrieved_var;
    }
    std::cout << "PREPAFM> The following parameters will be used for AFM setup:" << std::endl;
    std::cout << "******************************" << std::endl;
    for(auto it_afm=afm_opt.begin();it_afm!=afm_opt.end(); it_afm++) {
            std::cout << std::left 
                << std::setw(10) <<"[" + it_afm->first + "]" << " is set to " 
                << std::setw(10) << it_afm->second << std::endl;
    }
    std::cout << "******************************" << std::endl;
    /* afm initialization */
    fix_atom(afm_info.nterm);
    auto natom = get_psf_info().natom;
    afm_info.do_afm = true;
    afm_info.cterm = natom;
    afm_info.max_dist = (natom-2)*3.8;
    afm_info.v_afm = 1e-5;
    afm_info.k_afm = 200.0;
    auto i = afm_info.nterm - 1;
    auto j = afm_info.cterm - 1;
    Vector r_n(get_xcoor(i),get_ycoor(i),get_zcoor(i));
    Vector r_c(get_xcoor(j),get_ycoor(j),get_zcoor(j));
    Vector r_nc = r_n - r_c;
    afm_info.tmp_dist = r_nc.norm();
    /**********************/
    return true;
}

void apply_afm(real tstep) {
    auto i = afm_info.nterm - 1;
    auto j = afm_info.cterm - 1;
    Vector r_n(get_xcoor(i),get_ycoor(i),get_zcoor(i));
    Vector r_c(get_xcoor(j),get_ycoor(j),get_zcoor(j));
    Vector r_nc = r_n - r_c;
    real tmp_dist = r_nc.norm();
    if(tmp_dist > afm_info.max_dist) return;
    real k_afm = afm_info.k_afm;
    real v_afm = afm_info.v_afm;
    real force = k_afm * (tmp_dist - afm_info.tmp_dist);
    Vector f_i = -1.0 * force * r_nc.unitvec();
    Vector f_j = -1.0 * f_i;
    add_to_force(i, f_i);
    add_to_force(j, f_j);
    afm_info.tmp_dist += v_afm*tstep;
}

const Afm_Info& get_afm_info() {
    return afm_info;
}