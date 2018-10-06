#ifndef COOR_HPP
#define COOR_HPP

#include <vector>
#include "define.hpp"

class CorData {
    private:
        // data, for easier DCD writing.
        Int _sizeof_coor;
        std::vector<float> _xcoor;
        std::vector<float> _ycoor;
        std::vector<float> _zcoor;
    public:
        // setter
        bool read_cor(const Str& inp_name);
        void move_cor(Vec3d& dr, const Int& atomid, bool is_movable);
        // getter
        const Vec3d get_atom_coor(const Int& atomid) const;
        // getter of pointer for writing dcd
        const float* p_xcoor(void);
        const float* p_ycoor(void);
        const float* p_zcoor(void);
        void print(void) const;
        Int get_coorsize() const;
};

#endif
