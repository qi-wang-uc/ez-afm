#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <utility>
#include "../include/param.h"
#include "../include/util.h"

/* PRM DATA STORAGE */
static std::map<PrmBondType,PrmBondParam> bond_params;
static std::map<PrmAnglType,PrmAnglParam> angl_params;
static std::map<PrmNbndType,PrmNbndParam> nbnd_params;
static std::map<PrmNbfxType,PrmNbfxParam> nbfx_params;
/* dihedrals are treated differently because Fourier terms are in different lines */
static std::map<PrmDiheType,std::vector<PrmDiheParam> > dihe_params;
/********************/

bool read_prm(const std::string &inp_name) {
	std::ifstream inp_file(inp_name);
	if(!inp_file.is_open()) {
		std::cout << "ERROR> File not found!" << std::endl;
		return false;
	} else {
		std::cout << "READPRM> Reading energy parameters from file ["
			<< inp_name << "]" << std::endl;
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
	bool is_reading_bond = false;
	bool is_reading_angl = false;
	bool is_reading_dihe = false;
	bool is_reading_nbnd = false;
	bool is_reading_nbfx = false;
	/* Container for buffer streams of each energy term */
	std::vector<std::string> buffer_bond;
	std::vector<std::string> buffer_angl;
	std::vector<std::string> buffer_dihe;
	std::vector<std::string> buffer_nbnd;
	std::vector<std::string> buffer_nbfx;
	/* Now read the parameter file */
	std::string each_line;
	while(std::getline(inp_file, each_line)) {
		if(each_line.empty()) continue;
		if(each_line.substr(0,1).compare("*")==0) std::cout << each_line << std::endl;
		if(each_line.substr(0,1).compare("!")==0) continue;
		/* Title, empty or commented lines have already been skipped */
		if(is_ignore_case_equal(each_line.substr(0,4),"BOND")) {
			std::cout << "READPRM> Reading PARAM BOND" << std::endl;
			is_reading_bond = true;
			continue;
		}
		if(is_ignore_case_equal(each_line.substr(0,4),"ANGL") ||
		   is_ignore_case_equal(each_line.substr(0,4),"THET")) {
			std::cout << "READPRM> Reading PARAM ANGLE" << std::endl;
			is_reading_angl = true;
			is_reading_bond = false;
			continue;
		}
		if(is_ignore_case_equal(each_line.substr(0,4),"DIHE") ||
		   is_ignore_case_equal(each_line.substr(0,3),"PHI")) {
			std::cout << "READPRM> Reading PARAM DIHEDRAL" << std::endl;
			is_reading_dihe = true;
			is_reading_angl = false;
			continue;
		}
		if(is_ignore_case_equal(each_line.substr(0,4),"NBON")) {
			std::cout << "READPRM> Reading PARAM NONBOND" << std::endl;
			is_reading_nbnd = true;
			is_reading_dihe = false;
			continue;
		}
		if(is_ignore_case_equal(each_line.substr(0,4),"NBFI")) {
			std::cout << "READPRM> Reading PARAM NBFIX" << std::endl;
			is_reading_nbnd = false;
			is_reading_nbfx = true;
			continue;
		}
		if(is_ignore_case_equal(each_line.substr(0,3),"END")) {
			is_reading_nbfx = false;
			continue;
		}
		if(is_reading_bond) buffer_bond.push_back(each_line);
		if(is_reading_angl) buffer_angl.push_back(each_line);
		if(is_reading_dihe) buffer_dihe.push_back(each_line);
		if(is_reading_nbnd) buffer_nbnd.push_back(each_line);
		if(is_reading_nbfx) buffer_nbfx.push_back(each_line);
	}
	inp_file.close();
	// Processing buffers
	parse_PrmBond(buffer_bond);
	parse_PrmAngl(buffer_angl);
	parse_PrmDihe(buffer_dihe);
	parse_PrmNbnd(buffer_nbnd);
	parse_PrmNbfx(buffer_nbfx);
	return true;
}

void parse_PrmBond(const std::vector<std::string> &inp_data) {
	std::string i;
	std::string j;
	real k;
	real r;
	std::stringstream inp_stream;
	for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
		inp_stream.str(*it);
		inp_stream >> i >> j >> k >> r;
		bond_params.insert(std::make_pair(PrmBondType(i,j),PrmBondParam(k,r)));
		inp_stream.clear();
	}
}

void parse_PrmAngl(const std::vector<std::string> &inp_data) {
	std::string i;
	std::string j;
	std::string k;
	real k_angle;
	real theta_eq;
	std::stringstream inp_stream;
	for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
		inp_stream.str(*it);
		inp_stream >> i >> j >> k >> k_angle >> theta_eq;
		angl_params.insert(std::make_pair(PrmAnglType(i,j,k),PrmAnglParam(k_angle,theta_eq)));
		inp_stream.clear();
	}
}

void parse_PrmDihe(const std::vector<std::string> &inp_data) {
	std::string i;
	std::string j;
	std::string k;
	std::string l;
	real f;	// force constant
	real d;	// delta
	unsigned short m;	// multiplicity
	std::stringstream inp_stream;
	for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
		inp_stream.str(*it);
		inp_stream >> i >> j >> k >> l >> f >> m >> d;
		auto i_dihe_type = PrmDiheType(i,j,k,l);
		auto i_dihe_param = PrmDiheParam(f,m,d);
		if(dihe_params.find(i_dihe_type)!=dihe_params.end()) {
			dihe_params[i_dihe_type].push_back(i_dihe_param);
		} else {
			std::vector<PrmDiheParam> vec_dihe_param = {i_dihe_param};
			dihe_params[i_dihe_type] = vec_dihe_param;
		}
		inp_stream.clear();
	}
}

void parse_PrmNbnd(const std::vector<std::string> &inp_data) {
	std::string i;
	real d;
	real eps;
	real em;
	std::stringstream inp_stream;
	for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
		inp_stream.str(*it);
		inp_stream >> i >> d >> eps >> em;
		nbnd_params.insert(std::make_pair(PrmNbndType(i),PrmNbndParam(d,eps,em)));
		inp_stream.clear();
	}
}

void parse_PrmNbfx(const std::vector<std::string> &inp_data) {
	std::string i;
	std::string j;
	real e;
	real s;
	real c;
	real d;
	std::stringstream inp_stream;

	for(auto it=inp_data.cbegin();it!=inp_data.cend();++it) {
		inp_stream.str(*it);
		inp_stream >> i >> j >> e >> s >> c >> d;
		// The sign of coefficient C need to be changed because it sits in the
		// column of a CHARMM parameter file which is supposed to be negative.
		c = -1.0*c;
		nbfx_params.insert(std::make_pair(PrmNbfxType(i,j),PrmNbfxParam(e,s,c,d)));
		inp_stream.clear();
	}
}

PrmBondParam get_bond_params(PrmBondType type_query) {
	return bond_params[type_query];
};

PrmAnglParam get_angl_params(PrmAnglType type_query) {
	return angl_params[type_query];
};

const std::vector<PrmDiheParam>& get_dihe_params(PrmDiheType type_query) {
	return dihe_params[type_query];
};

PrmNbndParam get_nbnd_params(PrmNbndType type_query) {
	return nbnd_params[type_query];
};

PrmNbfxParam get_nbfx_params(PrmNbfxType type_query) {
	return nbfx_params[type_query];
};

bool is_good_params(){
	if(bond_params.empty())
		std::cout << "WARNING> No bond parameters (optional) found." << std::endl;
	if(bond_params.empty())
		std::cout << "WARNING> No angle parameters (optional) found." << std::endl;
	if(bond_params.empty())
		std::cout << "WARNING> No dihedral parameters (optional) found." << std::endl;
	if(bond_params.empty())
		std::cout << "WARNING> No nbonded parameters (optional) found." << std::endl;
	if(nbfx_params.empty()) {
		std::cout << "ERROR> No NBFIX parameters (required) found." << std::endl;
		return false;
	}
	return true;
}