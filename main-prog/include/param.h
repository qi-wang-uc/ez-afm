#ifndef PARAM_H
#define PARAM_H

#include <map>
#include <vector>
#include <string>
#include "../include/main.h"
/* 
	To keep consistent with CHARMM or Fortran coding style, here "BOND",
	"ANGL", "DIHE", "NBND" and "NBFX" is used for "bond", "angle", "dihedral",
	"non-bonded" and "non-bonded fix" terms.

	The "Type" structs are used as Keys of the std::map, thus one should define
	a comparison operator, in this case the "<" is overloaded.
*/

struct PrmBondType {
	std::string type_i;
	std::string type_j;
	PrmBondType() {}
	PrmBondType(std::string i, std::string j) {
		if(i <= j) {
			type_i = i;
			type_j = j;
		} else {
			type_i = j;
			type_j = i;
		}
	}
	bool operator < (const PrmBondType tmp) const {
		return (type_i+type_j) < (tmp.type_i+tmp.type_j);
	}
};

struct PrmBondParam {
	real k_bond;
	real r_eq;
	PrmBondParam() {}
	PrmBondParam(real k, real r) {
		k_bond = k;
		r_eq = r;
	}
};

struct PrmAnglType {
	std::string type_i;
	std::string type_j;
	std::string type_k;
	PrmAnglType() {}
	PrmAnglType(std::string i, std::string j, std::string k) {
		type_j = j;
		if(i <= k) {
			type_i = i;
			type_k = k;
		} else {
			type_i = k;
			type_k = i;
		}
	}
	bool operator < (const PrmAnglType tmp) const {
		return (type_i+type_j+type_k) < (tmp.type_i+tmp.type_j+tmp.type_k);
	}
};

struct PrmAnglParam {
	real k_angle;
	real theta_eq;
	PrmAnglParam() {}
	PrmAnglParam(real k, real t) {
		k_angle = k;
		theta_eq = t;
	}
};

struct PrmDiheType {
	std::string type_i;
	std::string type_j;
	std::string type_k;
	std::string type_l;
	PrmDiheType() {}
	PrmDiheType(std::string i, std::string j, std::string k, std::string l) {
		if(j <= k) {
			type_j = j;
			type_k = k;
		} else {
			type_j = k;
			type_k = j;
		}
		if(i <= l) {
			type_i = i;
			type_l = l;
		} else {
			type_i = l;
			type_l = i;
		}
	}
	bool operator < ( const PrmDiheType tmp) const {
		return (type_i+type_j+type_k+type_l) < (tmp.type_i+tmp.type_j+tmp.type_k+tmp.type_l);
	}
	bool operator == (const PrmDiheType tmp) const {
		return (type_i+type_j+type_k+type_l) == (tmp.type_i+tmp.type_j+tmp.type_k+tmp.type_l);
	}
};

struct PrmDiheParam {
	real k_dihe;
	unsigned short mul;
	real delta;
	PrmDiheParam() {}
	PrmDiheParam(real k, unsigned short m, real d) {
		k_dihe = k;
		mul = m;
		delta = d;
	}
	bool operator == (const PrmDiheParam tmp) const {
		return (k_dihe==tmp.k_dihe) && (mul==tmp.mul) && (delta==tmp.delta) ;
	}
};

struct PrmNbndType {  // This is ignored in current BLN model.
	std::string type_i;
	PrmNbndType() {}
	PrmNbndType(std::string t) {
		type_i = t;
	}
	bool operator < ( const PrmNbndType tmp) const {
		return type_i < tmp.type_i;
	}
};

struct PrmNbndParam {
	real dummy;	// ignored if eps in negative in CHARMM.
	real epsilon;		// well-depth of LJ potential.
	real emin;	// emin/2 in CHARMM.
	PrmNbndParam() {}
	PrmNbndParam(real d, real eps, real em) {
		dummy = d;
		epsilon = eps;
		emin = em;
	}
};

struct PrmNbfxType {
	std::string type_i;
	std::string type_j;
	PrmNbfxType() {}
	PrmNbfxType(std::string i, std::string j) {
		type_i = i;
		type_j = j;
	}
	bool operator < ( const PrmNbfxType tmp) const {
		return (type_i+type_j) < (tmp.type_i+tmp.type_j);
	}
};

struct PrmNbfxParam {
	real epsilon;
	real sigma;
	real coeff_c;  
	real coeff_d;
	PrmNbfxParam() {}
	PrmNbfxParam(real e, real s, real c, real d) {
		epsilon = e;
		sigma = s;
		coeff_c = c;
		coeff_d = d;
	}
};

bool read_prm(const std::string &inp_name);

/* For setting up parameters (setters) */
void parse_PrmBond(const std::vector<std::string> &inp_data);
void parse_PrmAngl(const std::vector<std::string> &inp_data);
void parse_PrmDihe(const std::vector<std::string> &inp_data);
void parse_PrmNbnd(const std::vector<std::string> &inp_data);
void parse_PrmNbfx(const std::vector<std::string> &inp_data);

/* For retrieving paramters (getters) */
PrmBondParam get_bond_params(PrmBondType type_query);
PrmAnglParam get_angl_params(PrmAnglType type_query);
const std::vector<PrmDiheParam>& get_dihe_params(PrmDiheType type_query);
PrmNbndParam get_nbnd_params(PrmNbndType type_query);
PrmNbfxParam get_nbfx_params(PrmNbfxType type_query);

/* For sanity check of energy calculation 
   in this case, only NBFIX is requred, other energy terms are optional.
   TODO: what is a good parameter file?
*/
bool is_good_params();

#endif
