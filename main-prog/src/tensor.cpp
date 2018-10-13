#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>

#include "../include/tensor.hpp"
#include "../include/util.hpp"

/* To be further investigated */
const Real hiRadius = 3.6;
const Real rSpheres = 3.6;

void HItensor::init(const Int& natom) {
    this->ld = 3*natom;
    const auto LD = this->ld;
    this->D.resize(LD); this->S.resize(LD);
    for(Int i=0; i<LD; ++i) {
        D[i].resize(LD); S[i].resize(LD);
    }
}

void HItensor::build(const CorData& cor, const Real& zeta) {
    for(auto& d: this->D) {
        std::fill(d.begin(), d.end(), 0.0);
    }
    const auto Natoms = cor.get_coorsize();
    const auto LD     = this->ld;
    const auto kTzeta = ROOMKBT/zeta;
    /* sigma_i,j should be variable in more sophisticated model.*/
    for (Int i=0; i<Natoms; ++i) {
        auto sigma_i = rSpheres;
    /* step1. Filling lower off-diagonal blocks of D. */
        for (Int j=0; j<i; ++j) {
            auto sigma_j = rSpheres;
            auto r_ij = cor.get_atom_coor(i) - cor.get_atom_coor(j);
            auto dist = r_ij.norm();
            auto fact = 1.0/r_ij.norm2();
            auto dmat = fact*(r_ij.tensor_product(r_ij));
            auto dij = (dist<2*hiRadius) ? calc_Dij_overlap(sigma_i, sigma_j, dist, kTzeta, dmat) 
                                         : calc_Dij_nonoverlap(sigma_i, sigma_j, dist, kTzeta, dmat);
            this->D[i*3+0][j*3+0] = dij.xx;
            this->D[i*3+0][j*3+1] = dij.xy;
            this->D[i*3+0][j*3+2] = dij.xz;
            this->D[i*3+1][j*3+0] = dij.yx;
            this->D[i*3+1][j*3+1] = dij.yy;
            this->D[i*3+1][j*3+2] = dij.yz;
            this->D[i*3+2][j*3+0] = dij.zx;
            this->D[i*3+2][j*3+1] = dij.zy;
            this->D[i*3+2][j*3+2] = dij.zz;
        }
    /* step2. Filling diagonal blocks of D. */
        this->D[i*3+0][i*3+0] = kTzeta;
        this->D[i*3+1][i*3+1] = kTzeta;
        this->D[i*3+2][i*3+2] = kTzeta;
    }
    /* step3. Filling upper off-diagonal blocks of D. */
    for(Int i=0; i<LD; ++i) {
        for(Int j=0; j<i; ++j)
            this->D[j][i] = this->D[i][j];
    }
}

void HItensor::cholesky() {
    for(auto& s: this->S) {
        std::fill(s.begin(), s.end(), 0.0);
    }
    const auto LD = this->ld;
    for (Int i=0; i<LD; ++i) {
        /* diagonal */
        Real s = 0.0;
        for (Int j=0; j<i; ++j) {
            s += this->D[j][i]*this->D[j][i];
        }
        this->S[i][i] = sqrt(this->D[i][i]-s);
        /* off-diagonal */
        for (Int k=i+1; k<LD; ++k) {
            s = 0.0;
            for (Int j=0; j<i; ++j) {
                s += this->D[j][k]*this->D[j][i];
            }
            this->S[i][k] = (this->D[i][k]-s)/this->D[i][i];
        }
    }
}

Mat3x3 calc_Dij_overlap(const Real& sigma_i, const Real& sigma_j, const Real& dist, const Real& kTzeta, const Mat3x3& dmat) {
    const auto sigma_sum = sigma_i + sigma_j;
    const auto sigma_avg = sigma_sum / 2;
    const auto alpha = 1.0 - 9.0 * dist / (32.0*sigma_avg);
    const auto beta  = 3.0 * dist / (32.0 * sigma_avg);
    return kTzeta*(alpha*I + beta*dmat);
}

Mat3x3 calc_Dij_nonoverlap(const Real& sigma_i, const Real& sigma_j, const Real& dist, const Real& kTzeta, const Mat3x3& dmat) {
    const auto sigma_sum = sigma_i + sigma_j;
    const auto sigma_avg = sigma_sum / 2;
    const auto sigma_sq  = sigma_avg*sigma_avg;
    const auto coeff     = 0.75*sigma_avg/dist;
    const auto alpha = coeff * (1.0 + 2.0*sigma_sq/(3.0*dist*dist));
    const auto beta  = coeff * (1.0 - 2.0*sigma_sq/(dist*dist));
    return kTzeta*(alpha*I + beta*dmat);
}

void HItensor::apply_disp_d(const Real& timecoeff, const PsfData& psf, const Energy& ener, CorData& cor) {
    const Int Natom = psf.get_NATOM();
    const Real coeffd = -1.0*timecoeff/ROOMKBT;
    
    for(Int i=0; i<Natom; ++i) {
        if(psf.is_movable(i)) {
            for (Int j=0; j<Natom; ++j) {
                cor.px()[i] += (this->D[i*3+0][j*3+0] * ener.px()[j] +
                                this->D[i*3+0][j*3+1] * ener.py()[j] +
                                this->D[i*3+0][j*3+2] * ener.pz()[j])*coeffd;
                cor.py()[i] += (this->D[i*3+1][j*3+0] * ener.px()[j] +
                                this->D[i*3+1][j*3+1] * ener.py()[j] +
                                this->D[i*3+1][j*3+2] * ener.pz()[j])*coeffd;
                cor.pz()[i] += (this->D[i*3+2][j*3+0] * ener.px()[j] +
                                this->D[i*3+2][j*3+1] * ener.py()[j] +
                                this->D[i*3+2][j*3+2] * ener.pz()[j])*coeffd;
            }
        }
    }
}

void HItensor::apply_disp_r(const Real& timecoeff, const PsfData& psf, const Energy& ener, CorData& cor, const Rand& rand) {
    const Int Natom = psf.get_NATOM();
    const Real coeffr = sqrt(2.0*timecoeff);
    
    for(Int i=0; i<Natom; ++i) {
        if(psf.is_movable(i)) {
            for (Int j=0; j<Natom; ++j) {
                cor.px()[i] += (this->S[i*3+0][j*3+0] * rand.px()[j] +
                                this->S[i*3+0][j*3+1] * rand.py()[j] +
                                this->S[i*3+0][j*3+2] * rand.pz()[j])*coeffr;
                cor.py()[i] += (this->S[i*3+1][j*3+0] * rand.px()[j] +
                                this->S[i*3+1][j*3+1] * rand.py()[j] +
                                this->S[i*3+1][j*3+2] * rand.pz()[j])*coeffr;
                cor.pz()[i] += (this->S[i*3+2][j*3+0] * rand.px()[j] +
                                this->S[i*3+2][j*3+1] * rand.py()[j] +
                                this->S[i*3+2][j*3+2] * rand.pz()[j])*coeffr;
            }
        }
    }
}

void HItensor::print_mat(const Str& mat_name) {
    std::cout << "Tensor> Begin" << std::endl;
    std::cout << std::fixed;
    if("D"==mat_name) {
        for(const auto& d : this->D) {
           for (const auto& e : d)
                std::cout << std::setprecision(5) << std::setw(10) << e << " ";
            std::cout << std::endl;
        }
    }
    if("S"==mat_name) {
        for(const auto& s : this->S) {
           for (const auto& e : s)
                std::cout << std::setprecision(5) << std::setw(10) << e << " ";
            std::cout << std::endl;
        }
    }
    std::cout << "Tensor> End" << std::endl << std::endl;
}

void HItensor::write_mat(const Str& mat_name) {
    std::ofstream out_file(mat_name);
    for(const auto& d : this->D) {
        for (const auto& e : d)
            out_file << std::fixed << std::setprecision(5) << std::setw(10) << e << " ";
        out_file << std::endl;
    }
    out_file.close();
}