#ifndef ENERGY_HPP
#define ENERGY_HPP

#include <vector>
#include <map>
#include <utility>
#include "define.hpp"
#include "psf.hpp"
#include "coor.hpp"
#include "param.hpp"
#include "afm.hpp"

using NonBondPair = std::pair<Int, Int>;

class Energy {
    private:
        RealVec _xgrad;
        RealVec _ygrad;
        RealVec _zgrad;
        std::map<NonBondPair, bool> _nonbond_table;
        // data
        Real _ebond     = 0.0;
        Real _eangle    = 0.0;
        Real _edihedral = 0.0;
        Real _enonbond  = 0.0;
        Real _etotal    = 0.0;
    public:
        // setter
        Int  update_nonbond(const Real& cutoff, const CorData& cor);
        void init_energy(const Int& N);
        void apply_force(const Int& atomid, const Vec3d& force);
        void compute_energy(const PsfData& psf, const PrmData& prm, const CorData& cor);
        Real ebond(const PsfData& psf, const PrmData& prm, const CorData& cor);
        Real eangle(const PsfData& psf, const PrmData& prm, const CorData& cor);
        Real edihedral(const PsfData& psf, const PrmData& prm, const CorData& cor);
        Real enonbond(const PsfData& psf, const PrmData& prm, const CorData& cor);
        template<typename T>
        void other_forces(const IntVec& coor_set, const T& force_set);
        // getter
        void   print_energy(const Int& istep, const bool& is_header) const;
        Vec3d  get_force(const Int& atomid) const;
};

template<typename T>
void Energy::other_forces(const IntVec& coor_set, const T& force_set) {
    if (coor_set.size()!=force_set.size()) return;
    auto it_coor  = coor_set.begin();
    auto it_force = force_set.begin();
    for(; it_coor!=coor_set.end()&&it_force!=force_set.end(); it_coor++, it_force++)
        this->apply_force(*it_coor, *it_force);
}

#endif