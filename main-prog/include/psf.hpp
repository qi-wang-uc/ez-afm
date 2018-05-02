#ifndef PSF_HPP
#define PSF_HPP

#include <string>
#include <vector>
#include "afm.hpp"

struct PsfAtom {
    /* need to be algined? */
    size_t atom_id;
    std::string seg_name;
    size_t resi_id;
    std::string resi_name;
    std::string atom_name;
    std::string atom_type;
    double charge;
    double mass;
    size_t imove;
    double e_negativity;
    double hardness;
};

struct PsfBond {
    size_t atom_i;
    size_t atom_j;
    PsfBond() {}
    PsfBond(size_t i, size_t j) {
        atom_i = i;
        atom_j = j;
    }
};

struct PsfAngle {
    size_t atom_i;
    size_t atom_j;
    size_t atom_k;
    PsfAngle() {}
    PsfAngle(size_t i, size_t j, size_t k) {
        atom_i = i;
        atom_j = j;
        atom_k = k;
    }
};

struct PsfDihedral {
    size_t atom_i;
    size_t atom_j;
    size_t atom_k;
    size_t atom_l;
    PsfDihedral() {}
    PsfDihedral(size_t i, size_t j, size_t k, size_t l) {
        atom_i = i;
        atom_j = j;
        atom_k = k;
        atom_l = l;
    }
};

/* CMAP and other topology information are currently not implemented*/
class PsfData {
    private:
        size_t _NATOM;
    protected:
        std::vector<PsfAtom> _psf_atom;
        std::vector<PsfBond> _psf_bond;
        std::vector<PsfAngle> _psf_angle;
        std::vector<PsfDihedral> _psf_dihedral;
    public:
    //setter
        bool read_psf(const std::string& inp_name);
        void set_PsfAtom(const std::string& inp_data, const size_t& max_lines);
        void set_PsfBond(const std::string& inp_data, const size_t& max_lines);
        void set_PsfAngle(const std::string& inp_data, const size_t& max_lines);
        void set_PsfDihedral(const std::string& inp_data, const size_t& max_lines);
        bool set_fix_atom (const size_t& index);
    //getter
        const std::string get_atom_type(const size_t& query_id) const;
        const bool is_movable (const size_t& index) const;
        const size_t get_NATOM() const;
};

// helper
const int NUM_ENTRY_BOND     = 4;
const int NUM_ENTRY_ANGLE    = 3;
const int NUM_ENTRY_DIHEDRAL = 2;
size_t get_psf_section_count(std::string& inp_str);

#endif