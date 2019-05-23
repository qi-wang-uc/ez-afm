#ifndef PARAM_HPP
#define PARAM_HPP

#include <vector>
#include <string>
#include <map>
#include "define.hpp"

struct PrmBondType {
    Str type_i;
    Str type_j;
    PrmBondType() {}
    PrmBondType(Str i, Str j):type_i(i<=j?i:j),type_j(i<=j?j:i) {}
    bool operator < (const PrmBondType tmp) const {
        return (type_i+type_j) < (tmp.type_i+tmp.type_j);
    }
};

struct PrmBondParam {
    Real k_bond = 0.0;
    Real r_eq   = 0.0;
    PrmBondParam() {}
    PrmBondParam(Real k, Real r):k_bond(k),r_eq(r) {}
};

struct PrmAngleType {
    Str type_i;
    Str type_j;
    Str type_k;
    PrmAngleType() {}
    PrmAngleType(Str i, Str j, Str k):type_i(i<=k?i:k),type_j(j),type_k(i<=k?k:i) {}
    bool operator < (const PrmAngleType tmp) const {
        return (type_i+type_j+type_k) < (tmp.type_i+tmp.type_j+tmp.type_k);
    }
};

struct PrmAngleParam {
    Real k_angle  = 0.0;
    Real theta_eq = 0.0;
    PrmAngleParam() {}
    PrmAngleParam(Real k, Real t):k_angle(k),theta_eq(t) {}
};

struct PrmDihedralType {
    Str type_i;
    Str type_j;
    Str type_k;
    Str type_l;
    PrmDihedralType() {}
    PrmDihedralType(Str i, Str j, Str k, Str l):type_i(i<=l?i:l),type_j(j<=k?j:k),type_k(j<=k?k:j),type_l(i<=l?l:i) {}
    bool operator < ( const PrmDihedralType tmp) const {
        return (type_i+type_j+type_k+type_l) < (tmp.type_i+tmp.type_j+tmp.type_k+tmp.type_l);
    }
    bool operator == (const PrmDihedralType tmp) const {
        return (type_i+type_j+type_k+type_l) == (tmp.type_i+tmp.type_j+tmp.type_k+tmp.type_l);
    }
};

struct PrmDihedralParam {
    Real k_dihe  = 0.0;
    Int  mul     = 0;
    Real delta   = 0.0;
    PrmDihedralParam() {}
    PrmDihedralParam(Real k, Int m, Real d):k_dihe(k),mul(m),delta(d) {}
    bool operator == (const PrmDihedralParam tmp) const {
        return (k_dihe==tmp.k_dihe) && (mul==tmp.mul) && (delta==tmp.delta) ;
    }
};

struct PrmVdwType {
    Str type_i;
    PrmVdwType() {}
    PrmVdwType(Str t):type_i(t) {}
    bool operator < ( const PrmVdwType tmp) const {
        return type_i < tmp.type_i;
    }
};

struct PrmVdwParam {
    Real dummy   = 0.0; // ignored if eps is negative in CHARMM.
    Real epsilon = 0.0; // well-depth of LJ potential.
    Real emin    = 0.0; // emin/2 in CHARMM.
    PrmVdwParam() {}
    PrmVdwParam(Real d, Real eps, Real em):dummy(d),epsilon(eps),emin(em) {}
};

struct PrmNbfixType {
    Str type_i;
    Str type_j;
    PrmNbfixType() {}
    PrmNbfixType(Str i, Str j):type_i(i),type_j(j) {}
    bool operator < ( const PrmNbfixType tmp) const {
        return (type_i+type_j) < (tmp.type_i+tmp.type_j);
    }
};

struct PrmNbfixParam {
    Real epsilon = 0.0;
    Real sigma   = 0.0;
    Real coeff_c = 0.0;  
    Real coeff_d = 0.0;
    PrmNbfixParam() {}
    PrmNbfixParam(Real e, Real s, Real c, Real d):epsilon(e),sigma(s),coeff_c(c),coeff_d(d) {}
};

class PrmData {
    private:
        // data
        std::map<PrmBondType, PrmBondParam>  _param_bond;
        std::map<PrmAngleType,PrmAngleParam> _param_angle;
        std::map<PrmVdwType,  PrmVdwParam>   _param_vdw;
        std::map<PrmNbfixType,PrmNbfixParam> _param_nbfix;
        std::map<PrmDihedralType,std::vector<PrmDihedralParam> > _param_dihedral;
    public:
        // setters
        bool read_prm(const Str &inp_name);
        void parse_PrmBond(const std::vector<Str> &inp_data);
        void parse_PrmAngle(const std::vector<Str> &inp_data);
        void parse_PrmDihedral(const std::vector<Str> &inp_data);
        void parse_PrmVdw(const std::vector<Str> &inp_data);
        void parse_PrmNbfix(const std::vector<Str> &inp_data);
        // getters
        PrmBondParam  get_bond_params(PrmBondType type_query) const;
        PrmAngleParam get_angle_params(PrmAngleType type_query) const;
        std::vector<PrmDihedralParam> get_dihedral_params(PrmDihedralType type_query) const;
        PrmVdwParam get_vdw_params(PrmVdwType type_query) const;
        bool is_nbfix_type(PrmNbfixType type_query) const;
        PrmNbfixParam get_nbfix_params(PrmNbfixType type_query) const;
};

#endif
