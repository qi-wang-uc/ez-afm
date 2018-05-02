#ifndef RAND_HPP
#define RAND_HPP

#include <vector>
#include "../include/define.hpp"

/*******************************************************************
               GAUSSIAN RANDOM NUMBER DATA 
*******************************************************************/

class Rand {
    protected:
        // data
        std::vector<double> _xrand;
        std::vector<double> _yrand;
        std::vector<double> _zrand;
        size_t _size;
    public:
        // setter
        void  init_rand(const size_t& natom);
        void  gen_rand(const double& mean, const double& dev);
        // getter
        Vec3d get_rand(const size_t& atomid) const;
        void  print_rand(void) const;
};

#endif