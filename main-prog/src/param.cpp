#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include "../include/param.hpp"
#include "../include/util.hpp"

/* PRM DATA STORAGE */
/********************/

bool PrmData::read_prm(const std::string &inp_name) {
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> File not found!" << std::endl;
        return false;
    } else {
        std::cout << "ReadPRM> Reading energy parameters from file: " << inp_name << std::endl;
    }
    /* 
    Flags to determine which container the buffer should go to.
    >Note< Here the assumption is made that the input PRM file is
    organized as follows by convention: 
        [BOND] 
        [ANGLE/THETA] 
        [DIHEDRAL/PHI]
        [NBONDED]
        [NBFIX]
    More sophiscated functions can be used if those energy terms are
    in random order.	
    */
    bool is_reading_bond  = false;
    bool is_reading_angle = false;
    bool is_reading_dihedral = false;
    bool is_reading_vdw   = false;
    bool is_reading_nbfix = false;
    /* Container for buffer streams of each energy term */
    StrVec buffer_bond;
    StrVec buffer_angle;
    StrVec buffer_dihedral;
    StrVec buffer_vdw;
    StrVec buffer_nbfix;
    /* Now read the parameter file */
    std::string each_line;
    while(std::getline(inp_file, each_line)) {
        if(each_line.empty()) continue;
        if(each_line.substr(0,1).compare("*")==0) std::cout << each_line << std::endl;
        if(each_line.substr(0,1).compare("!")==0) continue;
        /* Title, empty or commented lines have already been skipped */
        if(is_ignore_case_equal(each_line.substr(0,4),"BOND")) {
            std::cout << "ReadPRM> Reading PARAM BOND" << std::endl;
            is_reading_bond = true;
            continue;
        }
        if(is_ignore_case_equal(each_line.substr(0,4),"ANGL") ||
            is_ignore_case_equal(each_line.substr(0,4),"THET")) {
            std::cout << "ReadPRM> Reading PARAM ANGLE" << std::endl;
            is_reading_angle = true;
            is_reading_bond = false;
            continue;
        }
        if(is_ignore_case_equal(each_line.substr(0,4),"DIHE") ||
            is_ignore_case_equal(each_line.substr(0,3),"PHI")) {
            std::cout << "ReadPRM> Reading PARAM DIHEDRAL" << std::endl;
            is_reading_dihedral = true;
            is_reading_angle = false;
            continue;
        }
        if(is_ignore_case_equal(each_line.substr(0,4),"NBON")) {
            std::cout << "ReadPRM> Reading PARAM NONBOND" << std::endl;
            is_reading_vdw = true;
            is_reading_dihedral = false;
            continue;
        }
        if(is_ignore_case_equal(each_line.substr(0,4),"NBFI")) {
            std::cout << "ReadPRM> Reading PARAM NBFIX" << std::endl;
            is_reading_vdw = false;
            is_reading_nbfix = true;
            continue;
        }
        if(is_ignore_case_equal(each_line.substr(0,3),"END")) {
            is_reading_nbfix = false;
            continue;
        }
        if(is_reading_bond) buffer_bond.push_back(each_line);
        if(is_reading_angle) buffer_angle.push_back(each_line);
        if(is_reading_dihedral) buffer_dihedral.push_back(each_line);
        if(is_reading_vdw) buffer_vdw.push_back(each_line);
        if(is_reading_nbfix) buffer_nbfix.push_back(each_line);
    }
    inp_file.close();
    this->parse_PrmBond(buffer_bond);
    this->parse_PrmAngle(buffer_angle);
    this->parse_PrmDihedral(buffer_dihedral);
    this->parse_PrmVdw(buffer_vdw);
    this->parse_PrmNbfix(buffer_nbfix);
    return true;
}

void PrmData::parse_PrmBond(const StrVec& inp_data) {
    Str i;
    Str j;
    Real k;
    Real r;
    std::stringstream inp_stream;
    for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
        inp_stream.str(*it);
        inp_stream >> i >> j >> k >> r;
        this->_param_bond.insert(std::make_pair(PrmBondType(i,j),PrmBondParam(k,r)));
        inp_stream.clear();
    }
}

void PrmData::parse_PrmAngle(const StrVec& inp_data) {
    Str i;
    Str j;
    Str k;
    Real k_angle;
    Real theta_eq;
    std::stringstream inp_stream;
    for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
        inp_stream.str(*it);
        inp_stream >> i >> j >> k >> k_angle >> theta_eq;
        this->_param_angle.insert(std::make_pair(PrmAngleType(i,j,k),PrmAngleParam(k_angle,theta_eq)));
        inp_stream.clear();
    }
}

void PrmData::parse_PrmDihedral(const StrVec& inp_data) {
    Str i;
    Str j;
    Str k;
    Str l;
    Real f;	// force constant
    Real d;	// delta
    unsigned short m;	// multiplicity
    std::stringstream inp_stream;
    for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
        inp_stream.str(*it);
        inp_stream >> i >> j >> k >> l >> f >> m >> d;
        auto i_dihe_type = PrmDihedralType(i,j,k,l);
        auto i_dihe_param = PrmDihedralParam(f,m,d);
        if(this->_param_dihedral.find(i_dihe_type)!=this->_param_dihedral.end()) {
            this->_param_dihedral[i_dihe_type].push_back(i_dihe_param);
        } else {
            std::vector<PrmDihedralParam> vec_dihe_param = {i_dihe_param};
            this->_param_dihedral[i_dihe_type] = vec_dihe_param;
        }
        inp_stream.clear();
    }
}

void PrmData::parse_PrmVdw(const std::vector<std::string> &inp_data) {
    Str  i;
    Real d;
    Real eps;
    Real em;
    std::stringstream inp_stream;
    for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
        inp_stream.str(*it);
        inp_stream >> i >> d >> eps >> em;
        this->_param_vdw.insert(std::make_pair(PrmVdwType(i),PrmVdwParam(d,eps,em)));
        inp_stream.clear();
    }
}

void PrmData::parse_PrmNbfix(const std::vector<std::string> &inp_data) {
    Str  i;
    Str  j;
    Real e;
    Real s;
    Real c;
    Real d;
    std::stringstream inp_stream;

    for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
        inp_stream.str(*it);
        inp_stream >> i >> j >> e >> s >> c >> d;
        /* The sign of coefficient [c] needs to be changed because it sits in the
            column of a CHARMM parameter file which is supposed to be negative. */
        c *= -1.0;
        this->_param_nbfix.insert(std::make_pair(PrmNbfixType(i,j),PrmNbfixParam(e,s,c,d)));
        inp_stream.clear();
    }
}

PrmBondParam PrmData::get_bond_params(PrmBondType type_query) const {
    return this->_param_bond.at(type_query);
};

PrmAngleParam PrmData::get_angle_params(PrmAngleType type_query) const{
    return this->_param_angle.at(type_query);
};

std::vector<PrmDihedralParam> PrmData::get_dihedral_params(PrmDihedralType type_query) const {
    return this->_param_dihedral.at(type_query);
};

PrmVdwParam PrmData::get_vdw_params(PrmVdwType type_query) const {
    return this->_param_vdw.at(type_query);
};

PrmNbfixParam PrmData::get_nbfix_params(PrmNbfixType type_query) const {
    return this->_param_nbfix.at(type_query);
};

bool PrmData::is_nbfix_type(PrmNbfixType type_query) const {
    return this->_param_nbfix.find(type_query) != this->_param_nbfix.end();
}
