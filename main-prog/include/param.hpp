#ifndef PARAM_HPP
#define PARAM_HPP

#include <vector>
#include <string>
#include <map>

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
	double k_bond = 0.0;
	double r_eq   = 0.0;
	PrmBondParam() {}
	PrmBondParam(double k, double r) {
		k_bond = k;
		r_eq = r;
	}
};

struct PrmAngleType {
	std::string type_i;
	std::string type_j;
	std::string type_k;
	PrmAngleType() {}
	PrmAngleType(std::string i, std::string j, std::string k) {
		type_j = j;
		if(i <= k) {
			type_i = i;
			type_k = k;
		} else {
			type_i = k;
			type_k = i;
		}
	}
	bool operator < (const PrmAngleType tmp) const {
		return (type_i+type_j+type_k) < (tmp.type_i+tmp.type_j+tmp.type_k);
	}
};

struct PrmAngleParam {
	double k_angle  = 0.0;
	double theta_eq = 0.0;
	PrmAngleParam() {}
	PrmAngleParam(double k, double t) {
		k_angle = k;
		theta_eq = t;
	}
};

struct PrmDihedralType {
	std::string type_i;
	std::string type_j;
	std::string type_k;
	std::string type_l;
	PrmDihedralType() {}
	PrmDihedralType(std::string i, std::string j, std::string k, std::string l) {
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
	bool operator < ( const PrmDihedralType tmp) const {
		return (type_i+type_j+type_k+type_l) < (tmp.type_i+tmp.type_j+tmp.type_k+tmp.type_l);
	}
	bool operator == (const PrmDihedralType tmp) const {
		return (type_i+type_j+type_k+type_l) == (tmp.type_i+tmp.type_j+tmp.type_k+tmp.type_l);
	}
};

struct PrmDihedralParam {
	unsigned short mul = 0;
	double k_dihe  = 0.0;
	double delta   = 0.0;
	PrmDihedralParam() {}
	PrmDihedralParam(double k, unsigned short m, double d) {
		k_dihe = k;
		mul = m;
		delta = d;
	}
	bool operator == (const PrmDihedralParam tmp) const {
		return (k_dihe==tmp.k_dihe) && (mul==tmp.mul) && (delta==tmp.delta) ;
	}
};

struct PrmVdwType {
	std::string type_i;
	PrmVdwType() {}
	PrmVdwType(std::string t) {
		type_i = t;
	}
	bool operator < ( const PrmVdwType tmp) const {
		return type_i < tmp.type_i;
	}
};

struct PrmVdwParam {
	double dummy   = 0.0; // ignored if eps is negative in CHARMM.
	double epsilon = 0.0; // well-depth of LJ potential.
	double emin    = 0.0; // emin/2 in CHARMM.
	PrmVdwParam() {}
	PrmVdwParam(double d, double eps, double em) {
		dummy = d;
		epsilon = eps;
		emin = em;
	}
};

struct PrmNbfixType {
	std::string type_i;
	std::string type_j;
	PrmNbfixType() {}
	PrmNbfixType(std::string i, std::string j) {
		type_i = i;
		type_j = j;
	}
	bool operator < ( const PrmNbfixType tmp) const {
		return (type_i+type_j) < (tmp.type_i+tmp.type_j);
	}
};

struct PrmNbfixParam {
	double epsilon = 0.0;
	double sigma   = 0.0;
	double coeff_c = 0.0;  
	double coeff_d = 0.0;
	PrmNbfixParam() {}
	PrmNbfixParam(double e, double s, double c, double d) {
		epsilon = e;
		sigma = s;
		coeff_c = c;
		coeff_d = d;
	}
};


class PrmData {
	protected:
		// data
		std::map<PrmBondType,PrmBondParam> _param_bond;
		std::map<PrmAngleType,PrmAngleParam> _param_angle;
		std::map<PrmVdwType,PrmVdwParam> _param_vdw;
		std::map<PrmNbfixType,PrmNbfixParam> _param_nbfix;
		std::map<PrmDihedralType,std::vector<PrmDihedralParam> > _param_dihedral;
	public:
		// setters
		bool read_prm(const std::string &inp_name);
		void parse_PrmBond(const std::vector<std::string> &inp_data);
		void parse_PrmAngle(const std::vector<std::string> &inp_data);
		void parse_PrmDihedral(const std::vector<std::string> &inp_data);
		void parse_PrmVdw(const std::vector<std::string> &inp_data);
		void parse_PrmNbfix(const std::vector<std::string> &inp_data);
		// getters
		PrmBondParam get_bond_params(PrmBondType type_query) const;
		PrmAngleParam get_angle_params(PrmAngleType type_query) const;
		std::vector<PrmDihedralParam> get_dihedral_params(PrmDihedralType type_query) const;
		PrmVdwParam get_vdw_params(PrmVdwType type_query) const;
		bool is_nbfix_type(PrmNbfixType type_query) const;
		PrmNbfixParam get_nbfix_params(PrmNbfixType type_query) const;
};
#endif
