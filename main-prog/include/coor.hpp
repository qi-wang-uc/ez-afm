#ifndef COOR_HPP
#define COOR_HPP

#include <vector>
#include "define.hpp"

class CorData {
    protected:
        // data, for easier DCD writing.
        std::vector<float> _xcoor;
        std::vector<float> _ycoor;
        std::vector<float> _zcoor;
        size_t _sizeof_coor;
    public:
        // setter
        bool read_cor(const std::string& inp_name);
        void move_cor(Vec3d& dr, const size_t& atomid);
        // getter
        const Vec3d get_atom_coor(const size_t& atomid) const;
        // getter of pointer for writing dcd
        const float* p_xcoor(void);
        const float* p_ycoor(void);
        const float* p_zcoor(void);
        void print(void) const;
};

#endif
