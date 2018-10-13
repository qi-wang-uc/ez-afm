#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "../include/energy.hpp"
#include "../include/define.hpp"
#include "../include/util.hpp"

void Energy::print_energy(const Int& istep, const bool& is_header) const {
    if (is_header) {
        std::cout << "EnerTitle> " 
                  << std::setw(10) << "STEP" 
                  << std::setw(16) << "E_BOND" 
                  << std::setw(16) << "E_ANGLE" 
                  << std::setw(16) << "E_DIHEDRAL" 
                  << std::setw(16) << "E_NONBOND" 
                  << std::setw(16) << "E_TOTAL" 
                  << std::endl;
    } else {
        std::cout << "EnerInfo> " 
                  << std::setw(10) << istep 
                  << std::setw(16) << _ebond 
                  << std::setw(16) << _eangle 
                  << std::setw(16) << _edihedral
                  << std::setw(16) << _enonbond
                  << std::setw(16) << _etotal 
                  << std::endl;
    }
}

void Energy::init_energy(const Int& N) {
    if(!(_xgrad.empty() && _ygrad.empty() && _zgrad.empty())) return;
    _xgrad.resize(N);
    _ygrad.resize(N);
    _zgrad.resize(N);
}

void Energy::apply_force(const Int& atomid, const Vec3d& force) {
    _xgrad[atomid] += force.x;
    _ygrad[atomid] += force.y;
    _zgrad[atomid] += force.z;
}

void Energy::compute_energy(const PsfData& psf, const PrmData& prm, const CorData& cor) {
    /* zero out force array */
    std::fill(_xgrad.begin(),_xgrad.end(), 0.0);
    std::fill(_ygrad.begin(),_ygrad.end(), 0.0);
    std::fill(_zgrad.begin(),_zgrad.end(), 0.0);
    /* calculate each energy term */
    _ebond     = ebond(psf, prm, cor);
    _eangle    = eangle(psf, prm, cor);
    _edihedral = edihedral(psf, prm, cor);
    _enonbond  = enonbond(psf, prm, cor);
    _etotal = _ebond + _eangle + _edihedral + _enonbond;
}

Vec3d Energy::get_force(const Int& atomid) const {
    return Vec3d(_xgrad.at(atomid), _ygrad.at(atomid), _zgrad.at(atomid));
}

Real Energy::ebond(const PsfData& psf, const PrmData& prm, const CorData& cor) {
    Real e_bond = 0.0;
    for(const auto& each_bond : psf.get_bond()) {
        /* retrieve energy parameters */
        Int index_i = each_bond.atom_i-1;
        Int index_j = each_bond.atom_j-1;
        std::string  type_i = psf.get_atom_type(index_i);
        std::string  type_j = psf.get_atom_type(index_j);
        Real k_bond = prm.get_bond_params(PrmBondType(type_i,type_j)).k_bond;
        Real r_eq   = prm.get_bond_params(PrmBondType(type_i,type_j)).r_eq;
        /* retrieven bond component coordinates */
        auto r_i = cor.get_atom_coor(index_i);
        auto r_j = cor.get_atom_coor(index_j);
        auto r_ij = r_i - r_j;
        Real r_diff = r_ij.norm() - r_eq;
        /* calculate force and energy */
        Real force  = k_bond * r_diff;
        Real energy = force  * r_diff;
        auto f_i = -2.0f * force * r_ij.unitvec(); 
        auto f_j = -1.0f * f_i;
        /* add forces to gradient array */
        this->apply_force(index_i, f_i);
        this->apply_force(index_j, f_j);
        e_bond += energy;
    }
    return e_bond;
}

Real Energy::eangle(const PsfData& psf, const PrmData& prm, const CorData& cor) {
    Real e_angle = 0.0;
    for(const auto& each_angle : psf.get_angle()) {
        /* retrieve angle energy parameters */
        Int index_i = each_angle.atom_i-1;
        Int index_j = each_angle.atom_j-1;
        Int index_k = each_angle.atom_k-1;
        std::string  type_i = psf.get_atom_type(index_i);
        std::string  type_j = psf.get_atom_type(index_j);
        std::string  type_k = psf.get_atom_type(index_k);
        Real k_angle  = prm.get_angle_params(PrmAngleType(type_i,type_j,type_k)).k_angle;
        Real theta_eq = prm.get_angle_params(PrmAngleType(type_i,type_j,type_k)).theta_eq;
        theta_eq *= DEG2RAD;
        /* retrieven angle component coordinates */
        auto r_i = cor.get_atom_coor(index_i);
        auto r_j = cor.get_atom_coor(index_j);
        auto r_k = cor.get_atom_coor(index_k);
        auto r_ij = r_i - r_j;
        auto r_kj = r_k - r_j;
        Real cos_theta = r_ij.unitvec().dot_product(r_kj.unitvec());
        Real sin_theta = sqrt(1-cos_theta*cos_theta);
        Real theta = acos(cos_theta);
        /* calculate force and energy */
        Real force  = 0.1 * k_angle * (theta - theta_eq);
        Real energy = force * (theta - theta_eq);
        force /= sin_theta*0.5;
        auto f_i = (force/r_ij.norm())*(r_kj.unitvec()-cos_theta*r_ij.unitvec());
        auto f_k = (force/r_kj.norm())*(cos_theta*r_kj.unitvec()-r_ij.unitvec());
        auto f_j = -1.0f * (f_i + f_k);
        /* apply changes into force array */
        this->apply_force(index_i, f_i);
        this->apply_force(index_j, f_j);
        this->apply_force(index_k, f_k);
        e_angle += energy;
    }
    return e_angle;
}

Real Energy::edihedral(const PsfData& psf, const PrmData& prm, const CorData& cor) {
    Real e_dihe = 0.0;
    for(const auto& each_dihedral : psf.get_dihedral()) {
        /* retrieve atom information */
        Int index_i = each_dihedral.atom_i-1;
        Int index_j = each_dihedral.atom_j-1;
        Int index_k = each_dihedral.atom_k-1;
        Int index_l = each_dihedral.atom_l-1;
        std::string  type_i = psf.get_atom_type(index_i);
        std::string  type_j = psf.get_atom_type(index_j);
        std::string  type_k = psf.get_atom_type(index_k);
        std::string  type_l = psf.get_atom_type(index_l);
        auto r_i = cor.get_atom_coor(index_i);
        auto r_j = cor.get_atom_coor(index_j);
        auto r_k = cor.get_atom_coor(index_k);
        auto r_l = cor.get_atom_coor(index_l);
        auto r_ij = r_i - r_j; 
        auto r_jk = r_j - r_k; auto r_kj = r_k - r_j; 
        auto r_kl = r_k - r_l; auto r_lk = r_l - r_k;
        Real rjk2inv = 1.0f / r_jk.norm2();
        /* vectors R and S form the "eigen-angle" of dihedral such that cos(<R,S>)=cos(phi)*/
        auto R = r_ij - (r_ij.dot_product(r_kj.unitvec()))*r_kj.unitvec();
        auto S = r_lk - (r_lk.dot_product(r_kj.unitvec()))*r_kj.unitvec();
        Real cos_phi = R.unitvec().dot_product(S.unitvec());
        /* vectors M and N are used to determine the sign of phi using IUPAC convention */
        auto M = r_ij.cross_product(r_kj);
        auto N = r_kj.cross_product(r_kl);
        Real dp = r_kj.dot_product(M.cross_product(N));
        Real sign_of_phi = (dp > 0.0) ? 1.0 : -1.0;
        Real phi = sign_of_phi * acos(cos_phi);
        /* Since f_i+f_j+f_k+f_l=0, choose vector T such that f_j=-f_i+T and f_k=-f_l-T. */
        auto dihe_prm = prm.get_dihedral_params(PrmDihedralType(type_i,type_j,type_k,type_l));
        for(const auto& each_term : dihe_prm) {
            Real k_dihe = each_term.k_dihe ;
            Int mul    = each_term.mul;
            Real delta  = each_term.delta * DEG2RAD;
            Real energy = k_dihe * (1.0f + cos(mul*phi-delta));
            Real force  = -1.0f * k_dihe * mul * sin(mul*phi-delta);
            auto f_i = -1.0 * force * r_kj.norm()* M /M.norm2();
            auto f_l =  1.0 * force * r_kj.norm()* N /N.norm2();
            auto T   = rjk2inv*(r_ij.dot_product(r_kj)*f_i-(r_kl.dot_product(r_kj)*f_l));
            auto f_j = -1.0f * f_i + T;
            auto f_k = -1.0f * f_l - T;
            /* apply changes to force array */
            this->apply_force(index_i, f_i);
            this->apply_force(index_j, f_j);
            this->apply_force(index_k, f_k);
            this->apply_force(index_l, f_l);
            e_dihe += energy;
        }
    }
    return e_dihe;
}

Int Energy::update_nonbond(const Real& cutoff, const CorData& cor) {
    const Real RCUT2 = cutoff*cutoff;
    const Int NATOM = cor.get_coorsize();
    Int pair_list_counter = 0;
    for(Int index_i=0; index_i < NATOM-RESISEP; ++index_i) {
        auto r_i = cor.get_atom_coor(index_i);
        for(Int index_j=index_i+RESISEP; index_j < NATOM; ++index_j) {
            auto r_j = cor.get_atom_coor(index_j);
            auto r_ij = r_i - r_j;
            if (r_ij.norm2() < RCUT2) {
                this->_nonbond_table[std::make_pair(index_i,index_j)] = true;
                pair_list_counter++;
            }
        }
    }
    return pair_list_counter;
}

Real Energy::enonbond (const PsfData& psf, const PrmData& prm, const CorData& cor) {
    Real e_nonbond = 0.0;
    for(const auto& each_pair : this->_nonbond_table) {
        if (!each_pair.second) continue;
        auto index_i = each_pair.first.first;
        auto index_j = each_pair.first.second;
        auto r_i = cor.get_atom_coor(index_i);
        auto r_j = cor.get_atom_coor(index_j);
        auto r_ij = r_i - r_j;
        std::string type_i = psf.get_atom_type(index_i);
        std::string type_j = psf.get_atom_type(index_j);
        /* For VDW potential */
        auto nb_prm_i = prm.get_vdw_params(PrmVdwType(type_i));
        auto nb_prm_j = prm.get_vdw_params(PrmVdwType(type_j));
        Real epsilon = sqrt(nb_prm_i.epsilon * nb_prm_j.epsilon);
        Real sigma   = (nb_prm_i.emin + nb_prm_j.emin) * 0.5;
        Real sigma2  = sigma*sigma;
        Real fr2 = sigma2 / r_ij.norm2();
        Real fr6 = fr2 * fr2 * fr2;
        Real energy = 4.0*epsilon*fr6*(fr6-1.0);
        Real force = 48.0*epsilon*fr6*(fr6-0.5)/r_ij.norm2();
        auto f_i = force * r_ij;
        auto f_j = -1.0f * f_i;
        this->apply_force(index_i, f_i);
        this->apply_force(index_j, f_j);
        e_nonbond += energy;
        /* For NBFIX potential */
        auto query_type = PrmNbfixType(type_i, type_j);
        if (!prm.is_nbfix_type(query_type)) continue;
        auto nbfx_prm = prm.get_nbfix_params(query_type);
        epsilon = -1. * nbfx_prm.epsilon;
        sigma   = nbfx_prm.sigma;
        sigma2  = sigma*sigma;
        fr2 = sigma2 / r_ij.norm2();
        fr6 = fr2 * fr2 * fr2;
        energy = 4.0*epsilon*fr6*(fr6-1.0);
        force = 48.0*epsilon*fr6*(fr6-0.5)/r_ij.norm2();
        f_i = force * r_ij;
        f_j = -1.0f * f_i;
        this->apply_force(index_i, f_i);
        this->apply_force(index_j, f_j);
        e_nonbond += energy;
    } 
    return e_nonbond;
}

const Real* Energy::px() const {
    return this->_xgrad.data();
}
const Real* Energy::py() const {
    return this->_ygrad.data();
}
const Real* Energy::pz() const {
    return this->_zgrad.data();
}