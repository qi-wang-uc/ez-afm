#ifndef PSF_HPP
#define PSF_HPP

#include <string>
#include <vector>
#include "afm.hpp"

struct PsfAtom {
    /* need to be algined? */
    Int atom_id;
    Str seg_name;
    Int resi_id;
    Str resi_name;
    Str atom_name;
    Str atom_type;
    Real charge;
    Real mass;
    Int  imove;
    Real e_negativity;
    Real hardness;
};

struct PsfBond {
    Int atom_i;
    Int atom_j;
    PsfBond() {}
    PsfBond(Int i, Int j):atom_i(i), atom_j(j) {}
};

struct PsfAngle {
    Int atom_i;
    Int atom_j;
    Int atom_k;
    PsfAngle() {}
    PsfAngle(Int i, Int j, Int k):atom_i(i),atom_j(j),atom_k(k) {}
};

struct PsfDihedral {
    Int atom_i;
    Int atom_j;
    Int atom_k;
    Int atom_l;
    PsfDihedral() {}
    PsfDihedral(Int i, Int j, Int k, Int l):atom_i(i),atom_j(j),atom_k(k), atom_l(l) {}
};

/* CMAP and other topology information are currently not implemented*/
class PsfData {
    private:
        Int _NATOM;
        std::vector<PsfAtom> _psf_atom;
        std::vector<PsfBond> _psf_bond;
        std::vector<PsfAngle> _psf_angle;
        std::vector<PsfDihedral> _psf_dihedral;
    public:
    //setter
        bool read_psf(const Str& inp_name);
        void set_PsfAtom(const Str& inp_data, const Int& max_lines);
        void set_PsfBond(const Str& inp_data, const Int& max_lines);
        void set_PsfAngle(const Str& inp_data, const Int& max_lines);
        void set_PsfDihedral(const Str& inp_data, const Int& max_lines);
        bool set_fix_atom (const Int& index);
    //getter
        const Str get_atom_type(const Int& query_id) const;
        const bool is_movable (const Int& index) const;
        const Int get_NATOM() const;
        std::vector<PsfBond> const& get_bond() const;
        std::vector<PsfAngle> const& get_angle() const;
        std::vector<PsfDihedral> const& get_dihedral() const;
};

// helper
const Int NUM_ENTRY_BOND     = 4;
const Int NUM_ENTRY_ANGLE    = 3;
const Int NUM_ENTRY_DIHEDRAL = 2;
Int get_psf_section_count(Str& inp_str);

#endif