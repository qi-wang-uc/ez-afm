#ifndef DEFINE_HPP
#define DEFINE_HPP

#include <vector>
#include <string>
#include <cmath>
#include <map>

using Int  = size_t;
using Str  = std::string;
using Real = double;
using StrVec  = std::vector<std::string>;
using IntVec  = std::vector<Int>;
using RealVec = std::vector<Real>;
using StrTable = std::map<std::string,std::string>;
using RealMat = std::vector<std::vector<Real>>;

/********************************************/
const Int RESISEP  = 3;         // Residue separation to determine nonbonded pairs.
const Real DEG2RAD = 0.0174533; // Degree to radian
const Real KB      = 0.0019899; // Boltzman constant
const Real TIMFAC  = 0.0488882; // Convert AKMA time unit to picosecond.
const Real FCONV   = 69.478508; // Converts forces from (kcal/mol)/A to pN
const Real ROOMKBT = 0.5961573; // energy at room temperature (300K) in kcal/mol
/********************************************/

struct Mat3x3 {
    Real xx; Real xy; Real xz;
    Real yx; Real yy; Real yz;
    Real zx; Real zy; Real zz;

    Mat3x3() {}
    Mat3x3(Real k) :
        xx(k),xy(k),xz(k),yx(k),yy(k),yz(k),zx(k),zy(k),zz(k) {}
    Mat3x3(Real xx, Real xy, Real xz, Real yx, Real yy, Real yz, Real zx, Real zy, Real zz):
        xx(xx),xy(xy),xz(xz),yx(yx),yy(yy),yz(yz),zx(zx),zy(zy),zz(zz) {}

    /* scalar algebra */
    friend Mat3x3 operator * (const Real& fact, const Mat3x3& mat) {
        return Mat3x3(mat.xx*fact, mat.xy*fact, mat.xz*fact,
                      mat.yx*fact, mat.yy*fact, mat.yz*fact,
                      mat.zx*fact, mat.zy*fact, mat.zz*fact);
    }
    Mat3x3& operator *= (const Real &fact) {
        this->xx *= fact; this->xy *= fact; this->xz *= fact;
        this->yx *= fact; this->yy *= fact; this->yz *= fact;
        this->zx *= fact; this->zy *= fact; this->zz *= fact;
        return *this;
    }
    friend Mat3x3 operator / (const Mat3x3 &mat, const Real &fact) {
        return Mat3x3(mat.xx/fact, mat.xy/fact, mat.xz/fact,
                      mat.yx/fact, mat.yy/fact, mat.yz/fact,
                      mat.zx/fact, mat.zy/fact, mat.zz/fact);
    }
    Mat3x3& operator /= (const Real &fact) {
        this->xx /= fact; this->xy /= fact; this->xz /= fact;
        this->yx /= fact; this->yy /= fact; this->yz /= fact;
        this->zx /= fact; this->zy /= fact; this->zz /= fact;
        return *this;
    }
    /* Mat3x3 algebra */
    Mat3x3& operator += (const Mat3x3& rhs) {
        this->xx += rhs.xx; this->xy += rhs.xy; this->xz += rhs.xz;
        this->yx += rhs.yx; this->yy += rhs.yy; this->yz += rhs.yz;
        this->zx += rhs.zx; this->zy += rhs.zy; this->zz += rhs.zz;
        return *this;
    }
    Mat3x3& operator -= (const Mat3x3& rhs) {
        this->xx -= rhs.xx; this->xy -= rhs.xy; this->xz -= rhs.xz;
        this->yx -= rhs.yx; this->yy -= rhs.yy; this->yz -= rhs.yz;
        this->zx -= rhs.zx; this->zy -= rhs.zy; this->zz -= rhs.zz;
        return *this;
    }
    
    friend Mat3x3 operator + (const Mat3x3& mat1, const Mat3x3& mat2) {
        return Mat3x3(mat1.xx+mat2.xx, mat1.xy+mat2.xy, mat1.xz+mat2.xz,
                      mat1.yx+mat2.yx, mat1.yy+mat2.yy, mat1.yz+mat2.yz,
                      mat1.zx+mat2.zx, mat1.zy+mat2.zy, mat1.zz+mat2.zz);

    }
    friend Mat3x3 operator - (const Mat3x3& mat1, const Mat3x3& mat2) {
        return Mat3x3(mat1.xx-mat2.xx, mat1.xy-mat2.xy, mat1.xz-mat2.xz,
                      mat1.yx-mat2.yx, mat1.yy-mat2.yy, mat1.yz-mat2.yz,
                      mat1.zx-mat2.zx, mat1.zy-mat2.zy, mat1.zz-mat2.zz);

    }
    /* matrix multiplication is not used durnig simulation thus not implemented */
};

struct Vec3d {
    /* basic elements and constructor */
    Real x; Real y; Real z;
    Vec3d() {}
    Vec3d(Real x, Real y, Real z) : x(x), y(y), z(z) {}
    /* intrinsic properties */
    Real norm(void) const {
        return sqrt(x*x + y*y + z*z);
    }
    Real norm2(void) const {
        return x*x + y*y + z*z;
    }
    Vec3d unitvec(void) const {
        auto norm = sqrt(x*x + y*y + z*z);
        return Vec3d(x/norm,y/norm,z/norm);
    }
    /* scalar algebra */
    friend Vec3d operator * (const Real &fact, const Vec3d &vec) {
        return Vec3d(vec.x*fact, vec.y*fact, vec.z*fact);
    }
    Vec3d& operator *= (const Real &fact) {
        x *= fact;
        y *= fact;
        z *= fact;
        return *this;
    }
    friend Vec3d operator / (const Vec3d &vec, const Real &fact) {
        return Vec3d(vec.x/fact, vec.y/fact, vec.z/fact);
    }
    Vec3d& operator /= (const Real &fact) {
        x /= fact;
        y /= fact;
        z /= fact;
        return *this;
    }
    /* vector algebra */
    Vec3d& operator += (const Vec3d &vec) {
        x += vec.x;
        y += vec.y;
        z += vec.z;
        return *this;
    }
    Vec3d& operator -= (const Vec3d &vec) {
        x -= vec.x;
        y -= vec.y;
        z -= vec.z;
        return *this;
    }
    friend Vec3d operator + (const Vec3d &vec1, const Vec3d &vec2) {
        return Vec3d(vec1.x+vec2.x, vec1.y+vec2.y, vec1.z+vec2.z);
    }
    friend Vec3d operator - (const Vec3d &vec1, const Vec3d &vec2) {
        return Vec3d(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z);
    }
    Vec3d cross_product(const Vec3d &vec) {
        return Vec3d(y*vec.z-vec.y*z, vec.x*z-x*vec.z, x*vec.y-vec.x*y);
    }
    Real dot_product(const Vec3d &vec) {
        return x*vec.x + y*vec.y + z*vec.z;
    }
    Mat3x3 tensor_product(const Vec3d &vec) {
        return Mat3x3(x*vec.x, x*vec.y, x*vec.z,
                      y*vec.x, y*vec.y, y*vec.z,
                      z*vec.x, z*vec.y, z*vec.z);
    }
};

#endif