#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include "../include/psf.hpp"
#include "../include/util.hpp"

bool PsfData::read_psf(const Str& inp_name) {
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open file: "  << inp_name << std::endl;
        return false;
    } else {
        std::cout << "ReadPSF> Reading PSF information from: " << inp_name << std::endl;
    }
    Str each_line;
    std::stringstream each_stream;
	
    bool is_reading_atom = false;
    bool is_reading_bond = false;
    bool is_reading_angle = false;
    bool is_reading_dihedral = false;
    /* 
    Currently only the NATOM, NBOND, NTHETA and NPHI sections are utilized
    from PSF file, so it's helpful to set a reading limit to avoid unwanted info. 
    max_lines = CEIL(maPsfAtomx_members / members_per_line); where members_per_line for
    NATOM, NBOND, NTHETA(NANGLE) and NPHI(NDIHEDRAL) are 1, 4, 3, 2, respectively.
    */
    Int max_members;
    Int max_lines;
    while(std::getline(inp_file, each_line)) {
        each_stream.str(each_line);
        if(each_line.empty()) continue;
        if(0==each_line.substr(0,3).compare("PSF")) continue;
		if(0==each_line.substr(0,1).compare("*")) std::cout << each_line << std::endl;
        /* Turn on or off switches (boolean flags) */
        if(each_line.find("!NATOM")!=std::string::npos) {
            max_members = get_psf_section_count(each_line);
            std::cout << "ReadPSF> Reading PSF ATOM " << encap(max_members) << std::endl;
            max_lines = max_members;
            is_reading_atom = true;
            continue;
        }
        if(each_line.find("!NBOND")!=std::string::npos) {
            max_members = get_psf_section_count(each_line);
            std::cout << "ReadPSF> Reading PSF BOND " << encap(max_members) << std::endl;
            max_lines = max_members/4 + 1;
            is_reading_atom = false;
            is_reading_bond = true;
            continue;
        }
        if(each_line.find("!NTHETA")!=std::string::npos) {
            max_members = get_psf_section_count(each_line);
            std::cout << "ReadPSF> Reading PSF ANGLE " << encap(max_members) << std::endl;
            max_lines = max_members/3 + 1;
            is_reading_bond = false;
            is_reading_angle = true;
            continue;
        }
        if(each_line.find("!NPHI")!=std::string::npos) {
            max_members = get_psf_section_count(each_line);
            std::cout << "ReadPSF> Reading PSF DIHEDRAL " << encap(max_members) << std::endl;
            max_lines = max_members/2 + 1;
            is_reading_angle = false;
            is_reading_dihedral = true;
            continue;
        }
        // Now process each line based on flags.
        if(is_reading_atom) this->set_PsfAtom(each_line, max_lines);
        if(is_reading_bond) this->set_PsfBond(each_line, max_lines);
        if(is_reading_angle) this->set_PsfAngle(each_line, max_lines);
        if(is_reading_dihedral) this->set_PsfDihedral(each_line, max_lines);
    }
    inp_file.close();
    this->_NATOM = this->_psf_atom.size();
    return true;
}

void PsfData::set_PsfAtom(const Str& inp_data, const Int& max_lines) {
    static Int counter = 0;
    if(counter > max_lines) return;
    PsfAtom tmp_psf_atom;
    std::stringstream inp_stream(inp_data);
    inp_stream >> tmp_psf_atom.atom_id 
               >> tmp_psf_atom.seg_name 
               >> tmp_psf_atom.resi_id
               >> tmp_psf_atom.resi_name 
               >> tmp_psf_atom.atom_name 
               >> tmp_psf_atom.atom_type
               >> tmp_psf_atom.charge 
               >> tmp_psf_atom.mass 
               >> tmp_psf_atom.imove 
               >> tmp_psf_atom.e_negativity 
               >> tmp_psf_atom.hardness;
    this->_psf_atom.push_back(tmp_psf_atom);
    counter++;
}

void PsfData::set_PsfBond(const Str &inp_data, const Int& max_lines) {
    static Int counter = 0;
    if(counter > max_lines) return;
    std::stringstream inp_stream(inp_data);
    Int n_loop = n_of_words(inp_data)/2;
    n_loop = (n_loop==NUM_ENTRY_BOND) ? NUM_ENTRY_BOND : n_loop;
    for(Int i=0; i<n_loop; ++i) {
        PsfBond tmp_psf_bond;
        inp_stream >> tmp_psf_bond.atom_i 
                   >> tmp_psf_bond.atom_j;
        this->_psf_bond.push_back(tmp_psf_bond);
    }
    counter++;	
}

void PsfData::set_PsfAngle(const Str &inp_data, const Int& max_lines) {
    static Int counter = 0;
    if(counter > max_lines) return;
    std::stringstream inp_stream(inp_data);
    Int n_loop = n_of_words(inp_data)/3;
    n_loop = (n_loop==NUM_ENTRY_ANGLE) ? NUM_ENTRY_ANGLE : n_loop;
    for(Int i=0; i<n_loop; ++i) {
        PsfAngle tmp_psf_angle;
        inp_stream >> tmp_psf_angle.atom_i 
                   >> tmp_psf_angle.atom_j 
                   >> tmp_psf_angle.atom_k;
        this->_psf_angle.push_back(tmp_psf_angle);
    }
    counter++;
}

void PsfData::set_PsfDihedral(const Str &inp_data, const Int& max_lines) {
    static Int counter = 0;
    if(counter > max_lines) return;
    std::stringstream inp_stream(inp_data);
    Int n_loop = n_of_words(inp_data)/4;
    n_loop = (n_loop==NUM_ENTRY_DIHEDRAL) ? NUM_ENTRY_DIHEDRAL : n_loop;
    for(Int i=0; i<n_loop; ++i) {
        PsfDihedral tmp_psf_dihedral;
        inp_stream >> tmp_psf_dihedral.atom_i 
                   >> tmp_psf_dihedral.atom_j 
                   >> tmp_psf_dihedral.atom_k 
                   >> tmp_psf_dihedral.atom_l;
        this->_psf_dihedral.push_back(tmp_psf_dihedral);
    }
    counter++;
}

bool PsfData::set_fix_atom (const Int& index) {
    if(index > this->_psf_atom.size()) {
        std::cout << "ERROR> Requested fix-atom index " 
                  << encap(index)
                  << " exceeds maxium atom number."
                  << std::endl;
        return false;
	}
    std::cout << "FixAtom> Atom" << encap(index) << " will be fixed." << std::endl;
    this->_psf_atom[index-1].imove = 1;
    return true;
}

const Str PsfData::get_atom_type(const Int& query_id) const {
    return this->_psf_atom[query_id].atom_type;
}

const bool PsfData::is_movable (const Int& index) const {
    return 0 == this->_psf_atom[index].imove;
}

const Int PsfData::get_NATOM() const {
    return this->_NATOM;
}

Int get_psf_section_count(Str& inp_str) {
    auto npos = inp_str.find_first_not_of("0123456789");
    return std::stol(inp_str.substr(0,npos-1));
}

std::vector<PsfBond> const& PsfData::get_bond() const {
    return this->_psf_bond;
}

std::vector<PsfAngle> const& PsfData::get_angle() const {
    return this->_psf_angle;
}

std::vector<PsfDihedral> const& PsfData::get_dihedral() const {
    return this->_psf_dihedral;
}
