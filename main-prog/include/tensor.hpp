#ifndef TENSOR_HPP
#define TENSOR_HPP

#include "define.hpp"
#include "coor.hpp"
#include "psf.hpp"
#include "energy.hpp"
#include "rand.hpp"

/* A 3x3 REAL identity matrix */
const auto I = Mat3x3( 1.0, 0.0, 0.0, 
                       0.0, 1.0, 0.0, 
                       0.0, 0.0, 1.0 );

class HItensor {
    private:
        Int    ld; // leading dimension of D and S
        RealMat D; // Rotne-Prager-Yamakawa tensor
        RealMat S; // Tensor after cholesky decomposition
    public:
        HItensor() {}
        void init(const Int& Natom);
        void build(const CorData& cor, const Real& zeta);
        void gen_rand(const Real& mean, const Real& dev);
        void cholesky();
        void print_mat(const Str& mat_name);
        void write_mat(const Str& mat_name);
        // Determinant displacement and Random displacement
        void apply_disp_d(const Real& timecoeff, const PsfData& psf, const Energy& ener, CorData& cor);
        void apply_disp_r(const Real& timecoeff, const PsfData& psf, const Energy& ener, CorData& cor, const Rand& rand);
};

Mat3x3 calc_Dij_overlap(const Real& sigma_i, const Real& sigma_j, const Real& dist, const Real& kTzeta, const Mat3x3& dmat);
Mat3x3 calc_Dij_nonoverlap(const Real& sigma_i, const Real& sigma_j, const Real& dist, const Real& kTzeta, const Mat3x3& dmat);

#endif