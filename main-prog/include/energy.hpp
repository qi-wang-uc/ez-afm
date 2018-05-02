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

using NonBondPair = std::pair<size_t, size_t>;

class Energy : public PsfData, public CorData, public PrmData {
    protected:
        // data
        std::map<NonBondPair, bool> _nonbond_table;
        std::vector<double> _xgrad;
        std::vector<double> _ygrad;
        std::vector<double> _zgrad;
        double _ebond     = 0.0;
        double _eangle    = 0.0;
        double _edihedral = 0.0;
        double _enonbond  = 0.0;
        double _etotal    = 0.0;
    public:
        // setter
        void   init_energy(const size_t& N);
        size_t update_nonbond(const double& cutoff);
        void   apply_force(const size_t& atomid, const Vec3d& force);
        void   compute_energy(void);
        double ebond(void);
        double eangle(void);
        double edihedral(void);
        double enonbond(void);
        template<typename T>
        void other_forces(const std::vector<size_t>& coor_set, const T& force_set);
        // getter
        void   print_energy(const size_t& istep, const bool& is_header) const;
        Vec3d  get_force(const size_t& atomid) const;
};

template<typename T>
void Energy::other_forces(const std::vector<size_t>& coor_set, const T& force_set) {
    auto it_coor  = coor_set.begin();
    auto it_force = force_set.begin();
    if (coor_set.size()!=force_set.size()) return;
    for(; it_coor!=coor_set.end()&&it_force!=force_set.end(); it_coor++, it_force++)
        this->apply_force(*it_coor, *it_force);
}

#endif