#ifndef PSF_H
#define PSF_H

#include <string>
#include <vector>
#include "main.h"

/* For the ATOMs, urrently only the index and type are considered */
const int BOND_PER_LINE = 4;
const int ANGL_PER_LINE = 3;
const int DIHE_PER_LINE = 2;

struct PsfInfo {
    unsigned int natom;
    unsigned int nbond;
    unsigned int ntheta;
    unsigned int nphi;
};

struct PsfAtom {
    /* need to be algined? */
    unsigned int atom_id;
    std::string seg_name;
    unsigned int resi_id;
    std::string resi_name;
    std::string atom_name;
    std::string atom_type;
    real charge;
    real mass;
    unsigned short imove;
    real e_negativity;
    real hardness;
};

struct PsfBond {
    unsigned int atom_i;
    unsigned int atom_j;
    PsfBond() {}
    PsfBond(unsigned int i, unsigned int j) {
        atom_i = i;
        atom_j = j;
    }
};

struct PsfAngl {
    unsigned int atom_i;
    unsigned int atom_j;
    unsigned int atom_k;
    PsfAngl() {}
    PsfAngl(unsigned int i, unsigned int j, unsigned int k) {
        atom_i = i;
        atom_j = j;
        atom_k = k;
    }
};

struct PsfDihe {
    unsigned int atom_i;
    unsigned int atom_j;
    unsigned int atom_k;
    unsigned int atom_l;
    PsfDihe() {}
    PsfDihe(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
        atom_i = i;
        atom_j = j;
        atom_k = k;
        atom_l = l;
    }
};

/* CMAP and other topology information are currently not implemented*/

// main psf parsing function.
bool read_psf(const std::string &inp_name);

// parsing function of each section.
void parse_PsfAtom(const std::string &inp_data, const unsigned int& max_lines);
void parse_PsfBond(const std::string &inp_data, const unsigned int& max_lines);
void parse_PsfAngl(const std::string &inp_data, const unsigned int& max_lines);
void parse_PsfDihe(const std::string &inp_data, const unsigned int& max_lines);

// psf reading utility function.
unsigned int get_psf_section_count(std::string &inp_str);
const std::string get_atom_type(const unsigned int query_id);

bool fix_atom (unsigned int index);
bool is_movable (unsigned int index);

const std::vector<PsfBond>& psf_bond();
const std::vector<PsfAngl>& psf_angl();
const std::vector<PsfDihe>& psf_dihe();

// In this case, a good psf file is defined as non-empty file with at least 1 atom.
// TODO: what is a good psf file?
bool is_good_psf();
PsfInfo get_psf_info();

#endif
