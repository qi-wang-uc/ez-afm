#ifndef DEFINE_HPP
#define DEFINE_HPP

#include <cmath>

/********************************************/
const size_t RESISEP = 3;             // Residue separation to determine nonbonded pairs.
const double DEG2RAD = 0.0174533;     // Degree to radian
const double KB = 0.00198993;         // Boltzman constant
const double TIMFAC = 4.88882129e-02; // Convert AKMA time unit to picosecond.
const double FCONV  = 69.478508;      // Converts forces from (kcal/mol)/A to pN
/********************************************/

struct Vec3d {
    /* basic elements and constructor */
    double x; double y; double z;
    Vec3d() {}
    Vec3d(double x, double y, double z) : x(x), y(y), z(z) {}
    /* intrinsic properties */
    double norm(void) const {
        return sqrt(x*x + y*y + z*z);
    }
    double norm2(void) const {
        return x*x + y*y + z*z;
    }
    Vec3d unitvec(void) const {
        return Vec3d(x,y,z)/norm();
    }
    /* scalar algebra */
	friend Vec3d operator * (const double &fact, const Vec3d &vec) {
        return Vec3d(vec.x*fact, vec.y*fact, vec.z*fact);
    }
    Vec3d& operator *= (const double &fact) {
        x *= fact;
        y *= fact;
        z *= fact;
        return *this;
    }
	friend Vec3d operator / (const Vec3d &vec, const double &fact) {
        return Vec3d(vec.x/fact, vec.y/fact, vec.z/fact);
    }
    Vec3d& operator /= (const double &fact) {
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
    double dot_product(const Vec3d &vec) {
        return x*vec.x + y*vec.y + z*vec.z;
    }
};

#endif