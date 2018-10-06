#ifndef RAND_HPP
#define RAND_HPP

#include <vector>
#include "define.hpp"

/*******************************************************************
               GAUSSIAN RANDOM NUMBER DATA 
*******************************************************************/

class Rand {
    private:
        // data
        Int     _size;
        RealVec _xrand;
        RealVec _yrand;
        RealVec _zrand;
    public:
        // setter
        void  init_rand(const Int& natom);
        void  gen_rand(const Real& mean, const Real& dev);
        // getter
        Vec3d get_rand(const Int& atomid) const;
        void  print_rand(void) const;
};

#endif