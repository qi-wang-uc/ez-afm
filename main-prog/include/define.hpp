#ifndef DEFINE_HPP
#define DEFINE_HPP

#include <vector>
#include <string>
#include <cmath>
#include <map>

using Int  = int64_t;
using Str  = std::string;
using Real = double;
using StrVec  = std::vector<std::string>;
using IntVec  = std::vector<Int>;
using RealVec = std::vector<Real>;
using StrTable = std::map<std::string,std::string>;

/********************************************/
const Int RESISEP  = 3;              // Residue separation to determine nonbonded pairs.
const Real DEG2RAD = 0.0174533;      // Degree to radian
const Real KB      = 0.00198993;     // Boltzman constant
const Real TIMFAC  = 4.88882129e-02; // Convert AKMA time unit to picosecond.
const Real FCONV   = 69.478508;      // Converts forces from (kcal/mol)/A to pN
/********************************************/

struct Vec3d {
    /* basic elements and constructor */
    Real x;
    Real y;
    Real z;
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
    /* Vec3d algebra */
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
};

#endif