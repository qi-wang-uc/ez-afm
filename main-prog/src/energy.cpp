#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "../include/energy.hpp"
#include "../include/define.hpp"
#include "../include/util.hpp"

void Energy::print_energy(const size_t& istep, const bool& is_header) const {
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

void Energy::init_energy(const size_t& N) {
    if(!(_xgrad.empty() && _ygrad.empty() && _zgrad.empty())) return;
    _xgrad.resize(N);
    _ygrad.resize(N);
    _zgrad.resize(N);
}

void Energy::apply_force(const size_t& atomid, const Vec3d& force) {
    _xgrad[atomid] += force.x;
    _ygrad[atomid] += force.y;
    _zgrad[atomid] += force.z;
}

void Energy::compute_energy(void) {
    /* zero out force array */
    std::fill(_xgrad.begin(),_xgrad.end(), 0.0);
    std::fill(_ygrad.begin(),_ygrad.end(), 0.0);
    std::fill(_zgrad.begin(),_zgrad.end(), 0.0);
    /* calculate each energy term */
    _ebond     = ebond();
    _eangle    = eangle();
    _edihedral = edihedral();
    _enonbond  = enonbond();
    _etotal = _ebond + _eangle + _edihedral + _enonbond;
}

Vec3d Energy::get_force(const size_t& atomid) const {
    return Vec3d(_xgrad.at(atomid), _ygrad.at(atomid), _zgrad.at(atomid));
}

double Energy::ebond(void) {
    double e_bond = 0.0;
    for(const auto& each_bond : this->_psf_bond) {
        /* retrieve energy parameters */
        size_t index_i = each_bond.atom_i-1;
        size_t index_j = each_bond.atom_j-1;
        std::string  type_i = get_atom_type(index_i);
        std::string  type_j = get_atom_type(index_j);
        double k_bond = get_bond_params(PrmBondType(type_i,type_j)).k_bond;
        double r_eq   = get_bond_params(PrmBondType(type_i,type_j)).r_eq;
        /* retrieven bond component coordinates */
        auto r_i = get_atom_coor(index_i);
        auto r_j = get_atom_coor(index_j);
        auto r_ij = r_i - r_j;
        double r_diff = r_ij.norm() - r_eq;
        /* calculate force and energy */
        double force  = k_bond * r_diff;
        double energy = force  * r_diff;
        auto f_i = -2.0f * force * r_ij.unitvec(); 
        auto f_j = -1.0f * f_i;
        /* add forces to gradient array */
        this->apply_force(index_i, f_i);
        this->apply_force(index_j, f_j);
        e_bond += energy;
    }
    return e_bond;
}

double Energy::eangle (void) {
    double e_angle = 0.0;
    for(const auto& each_angle : this->_psf_angle) {
        /* retrieve angle energy parameters */
        size_t index_i = each_angle.atom_i-1;
        size_t index_j = each_angle.atom_j-1;
        size_t index_k = each_angle.atom_k-1;
        std::string  type_i = get_atom_type(index_i);
        std::string  type_j = get_atom_type(index_j);
        std::string  type_k = get_atom_type(index_k);
        double k_angle  = get_angle_params(PrmAngleType(type_i,type_j,type_k)).k_angle;
        double theta_eq = get_angle_params(PrmAngleType(type_i,type_j,type_k)).theta_eq;
        theta_eq *= DEG2RAD;
        /* retrieven angle component coordinates */
        auto r_i = get_atom_coor(index_i);
        auto r_j = get_atom_coor(index_j);
        auto r_k = get_atom_coor(index_k);
        auto r_ij = r_i - r_j;
        auto r_kj = r_k - r_j;
        double cos_theta = r_ij.unitvec().dot_product(r_kj.unitvec());
        double sin_theta = sqrt(1-cos_theta*cos_theta);
        double theta = acos(cos_theta);
        /* calculate force and energy */
        double force  = 0.1 * k_angle * (theta - theta_eq);
        double energy = force * (theta - theta_eq);
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

double Energy::edihedral(void) {
    double e_dihe = 0.0;
    for(const auto& each_dihedral : this->_psf_dihedral) {
        /* retrieve atom information */
        size_t index_i = each_dihedral.atom_i-1;
        size_t index_j = each_dihedral.atom_j-1;
        size_t index_k = each_dihedral.atom_k-1;
        size_t index_l = each_dihedral.atom_l-1;
        std::string  type_i = get_atom_type(index_i);
        std::string  type_j = get_atom_type(index_j);
        std::string  type_k = get_atom_type(index_k);
        std::string  type_l = get_atom_type(index_l);
        auto r_i = get_atom_coor(index_i);
        auto r_j = get_atom_coor(index_j);
        auto r_k = get_atom_coor(index_k);
        auto r_l = get_atom_coor(index_l);
        auto r_ij = r_i - r_j; 
        auto r_jk = r_j - r_k; auto r_kj = r_k - r_j; 
        auto r_kl = r_k - r_l; auto r_lk = r_l - r_k;
        double rjk2inv = 1.0f / r_jk.norm2();
        /* vectors R and S form the "eigen-angle" of dihedral such that cos(<R,S>)=cos(phi)*/
        auto R = r_ij - (r_ij.dot_product(r_kj.unitvec()))*r_kj.unitvec();
        auto S = r_lk - (r_lk.dot_product(r_kj.unitvec()))*r_kj.unitvec();
        double cos_phi = R.unitvec().dot_product(S.unitvec());
        /* vectors M and N are used to determine the sign of phi using IUPAC convention */
        auto M = r_ij.cross_product(r_kj);
        auto N = r_kj.cross_product(r_kl);
        double dp = r_kj.dot_product(M.cross_product(N));
        double sign_of_phi = (dp > 0.0) ? 1.0 : -1.0;
        double phi = sign_of_phi * acos(cos_phi);
        /* Since f_i+f_j+f_k+f_l=0, choose vector T such that f_j=-f_i+T and f_k=-f_l-T. */
        auto dihe_prm = get_dihedral_params(PrmDihedralType(type_i,type_j,type_k,type_l));
        for(const auto& each_term : dihe_prm) {
            double k_dihe = each_term.k_dihe ;
            size_t mul    = each_term.mul;
            double delta  = each_term.delta * DEG2RAD;
            double energy = k_dihe * (1.0f + cos(mul*phi-delta));
            double force  = -1.0f * k_dihe * mul * sin(mul*phi-delta);
            Vec3d f_i = -1.0 * force * r_kj.norm()* M /M.norm2();
            Vec3d f_l =  1.0 * force * r_kj.norm()* N /N.norm2();
            Vec3d T   = rjk2inv*(r_ij.dot_product(r_kj)*f_i-(r_kl.dot_product(r_kj)*f_l));
            Vec3d f_j = -1.0f * f_i + T;
            Vec3d f_k = -1.0f * f_l - T;
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

size_t Energy::update_nonbond(const double& cutoff) {
    const double RCUT2 = cutoff*cutoff;
    const size_t NATOM = this->get_NATOM();
    size_t pair_list_counter = 0;
    for(size_t index_i=0; index_i < NATOM-RESISEP; ++index_i) {
        auto r_i = get_atom_coor(index_i);
        for(size_t index_j=index_i+RESISEP; index_j < NATOM; ++index_j) {
            auto r_j = get_atom_coor(index_j);
            auto r_ij = r_i - r_j;
            if (r_ij.norm2() < RCUT2) {
                this->_nonbond_table[std::make_pair(index_i,index_j)] = true;
                pair_list_counter++;
            }
        }
    }
    return pair_list_counter;
}

double Energy::enonbond (void) {
    double e_nonbond = 0.0;
    for(const auto& each_pair : this->_nonbond_table) {
        if (!each_pair.second) continue;
        auto index_i = each_pair.first.first;
        auto index_j = each_pair.first.second;
        auto r_i = get_atom_coor(index_i);
        auto r_j = get_atom_coor(index_j);
        auto r_ij = r_i - r_j;
        std::string type_i = get_atom_type(index_i);
        std::string type_j = get_atom_type(index_j);
        /* For VDW potential */
        auto nb_prm_i = get_vdw_params(PrmVdwType(type_i));
        auto nb_prm_j = get_vdw_params(PrmVdwType(type_j));
        double epsilon = sqrt(nb_prm_i.epsilon * nb_prm_j.epsilon);
        double sigma   = (nb_prm_i.emin + nb_prm_j.emin) * 0.5;
        double sigma2  = sigma*sigma;
        double fr2 = sigma2 / r_ij.norm2();
        double fr6 = fr2 * fr2 * fr2;
        double energy = 4.0*epsilon*fr6*(fr6-1.0);
        double force = 48.0*epsilon*fr6*(fr6-0.5)/r_ij.norm2();
        auto f_i = force * r_ij;
        auto f_j = -1.0f * f_i;
        this->apply_force(index_i, f_i);
        this->apply_force(index_j, f_j);
        e_nonbond += energy;
        /* For NBFIX potential */
        auto query_type = PrmNbfixType(type_i, type_j);
        if (!is_nbfix_type(query_type)) continue;
        auto nbfx_prm = get_nbfix_params(query_type);
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

