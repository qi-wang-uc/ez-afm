#ifndef MAIN_H
#define MAIN_H

#include <cmath>

/* Perhaps it's not a good idea to define metrics here... */
using real = float;

struct Vector {
    /* basic elements and constructor */
    real x; real y; real z;
    Vector() {}
    Vector(real x, real y, real z) : x(x), y(y), z(z) {}
    /* intrinsic properties */
    real norm(void) const {
        return sqrt(x*x + y*y + z*z);
    }
    real norm2(void) const {
        return x*x + y*y + z*z;
    }
    Vector unitvec(void) const {
        return Vector(x,y,z)/norm();
    }
    /* scalar algebra */
	friend Vector operator*(const real &fact, const Vector &vec) {
        return Vector(vec.x*fact, vec.y*fact, vec.z*fact);
    }
    Vector& operator*=(const real &fact) {
        x *= fact;
        y *= fact;
        z *= fact;
        return *this;
    }
	friend Vector operator/(const Vector &vec, const real &fact) {
        return Vector(vec.x/fact, vec.y/fact, vec.z/fact);
    }
    Vector& operator/=(const real &fact) {
        x /= fact;
        y /= fact;
        z /= fact;
        return *this;
    }
    /* vector algebra */
    Vector& operator+= (const Vector &vec) {
        x += vec.x;
        y += vec.y;
        z += vec.z;
        return *this;
    }
    Vector& operator-= (const Vector &vec) {
        x -= vec.x;
        y -= vec.y;
        z -= vec.z;
        return *this;
    }
	friend Vector operator+(const Vector &vec1, const Vector &vec2) {
        return Vector(vec1.x+vec2.x, vec1.y+vec2.y, vec1.z+vec2.z);
    }
	friend Vector operator-(const Vector &vec1, const Vector &vec2) {
        return Vector(vec1.x-vec2.x, vec1.y-vec2.y, vec1.z-vec2.z);
    }
    Vector cross_product(const Vector &vec) {
        return Vector(y*vec.z-vec.y*z, vec.x*z-x*vec.z, x*vec.y-vec.x*y);
    }
    real dot_product(const Vector &vec) {
        return x*vec.x + y*vec.y + z*vec.z;
    }
};

struct Tensor {
    real xx; real xy; real xz;
    real yx; real yy; real yz;
    real zx; real zy; real zz;
};

#endif