#include <iostream>
#include <iomanip>
#include <map>
#include "../include/define.hpp"
#include "../include/energy.hpp"
#include "../include/coor.hpp"
#include "../include/uservar.hpp"

void AFM::setup_afm(StrVec cmds, UserVar& user_var) {
    std::cout << "PrepAFM> Preparing for AFM-like simulation setup." << std::endl;
    StrTable afm_opt;
    for(auto it_cmd=cmds.begin()+1; it_cmd!=cmds.end(); it_cmd+=2) {
        auto retrieved_var = user_var.query(*(it_cmd+1));
        afm_opt[*it_cmd] = retrieved_var;
    }
    std::cout << "PrepAFM> The following parameters will be used:" << std::endl;
    for(auto it_afm=afm_opt.begin();it_afm!=afm_opt.end(); it_afm++) {
        std::cout << std::left 
              << std::setw(10) << it_afm->first  << " = " 
              << std::setw(10) << it_afm->second
              << std::endl;
    }
    /* afm initialization */
    this->_afm_config.do_afm = true;
    this->_afm_config.nterm  = std::stol(afm_opt["nterm"]);
    this->_afm_config.cterm  = std::stol(afm_opt["cterm"]);
    this->_afm_config.k_afm  = std::stof(afm_opt["force"]);
    this->_afm_config.v_afm  = std::stof(afm_opt["velocity"]);
    this->_afm_config.max_dist = std::stof(afm_opt["maxdist"]);
    std::cout << "AfmTitle> "
              << std::setw(16) << "Force"
              << std::setw(16) << "Distance" 
              << std::endl;
}

AfmPair AFM::apply_afm(const Real& tstep, const AfmPair& afm_coors) {
    Vec3d r_nc = afm_coors[1] - afm_coors[0];
    Real tmp_dist = r_nc.norm();
    if(tmp_dist > _afm_config.max_dist)
        return AfmPair{Vec3d(0, 0, 0), Vec3d(0, 0, 0)};
    Real k_afm = _afm_config.k_afm;
    Real v_afm = _afm_config.v_afm;
    Real force = k_afm * (tmp_dist - _afm_config.tmp_dist);
    Vec3d f_i = -1.0 * force * r_nc.unitvec();
    Vec3d f_j = -1.0 * f_i;
    _afm_config.tmp_dist += v_afm*tstep;
    std::cout << "AfmInfo> " << std::fixed 
              << std::setw(16) << std::setprecision(6) << f_i.norm() 
              << std::setprecision(6) << std::setw(16) << tmp_dist 
              << std::endl;
    return AfmPair {f_i, f_j};
}

AfmConfig const& AFM::get_config() const {
    return this->_afm_config;
}