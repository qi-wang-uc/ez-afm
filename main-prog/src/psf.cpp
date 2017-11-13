#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include "../include/psf.h"
#include "../include/util.h"

static PsfInfo psf_info;
/* ATOM section are included in PSF file, but in this program it serves as atom-type
   and atom-movability query hash table. */
static std::map<unsigned int, PsfAtom> atom_topology;
/* PSF DATA STORAGE */
static std::vector<PsfBond> bond_topology;
static std::vector<PsfAngl> angl_topology;
static std::vector<PsfDihe> dihe_topology;
/********************/

bool read_psf(const std::string &inp_name) {
	std::ifstream inp_file(inp_name);
	if(!inp_file.is_open()) {
		std::cout << "ERROR> Cannot open file [" << inp_name << "]" << std::endl;
		return false;
	} else {
		std::cout << "READPSF> Reading PSF information from [" << inp_name 
			<< "]" << std::endl;
	}
	std::string each_line;
	std::stringstream each_stream;
	bool is_reading_atom = false;
	bool is_reading_bond = false;
	bool is_reading_angl = false;
	bool is_reading_dihe = false;
	/* 
	Currently only the NATOM, NBOND, NTHETA and NPHI sections are utilized
	from PSF file, so it's helpful to set a reading limit to avoid unwanted info. 
	max_lines = CEIL(maPsfAtomx_members / members_per_line); where members_per_line for
	NATOM, NBOND, NTHETA(NANGLE) and NPHI(NDIHEDRAL) are 1, 4, 3, 2, respectively.
	*/
	unsigned int max_members;
	unsigned int max_lines;
	while(std::getline(inp_file, each_line)) {
		each_stream.str(each_line);
		if(each_line.empty()) continue;
		if(each_line.substr(0,3).compare("PSF")==0) continue;
		if(each_line.substr(0,1).compare("*")==0) std::cout << each_line << std::endl;
		/* Turn on or off switches (boolean flags) */
		if(each_line.find("!NATOM")!=std::string::npos) {
			max_members = get_psf_section_count(each_line);
			psf_info.natom = max_members;	// Used for sanity check.
			std::cout << "READPSF> Reading PSF ATOM (" << max_members << ")" << std::endl;
			max_lines = max_members;
			is_reading_atom = true;
			continue;
		}
		if(each_line.find("!NBOND")!=std::string::npos) {
			max_members = get_psf_section_count(each_line);
			psf_info.nbond = max_members;
			std::cout << "READPSF> Reading PSF BOND (" << max_members << ")" << std::endl;
			max_lines = max_members/4 + 1;
			is_reading_atom = false;
			is_reading_bond = true;
			continue;
		}
		if(each_line.find("!NTHETA")!=std::string::npos) {
			max_members = get_psf_section_count(each_line);
			psf_info.ntheta = max_members;
			std::cout << "READPSF> Reading PSF ANGLE (" << max_members << ")" << std::endl;
			max_lines = max_members/3 + 1;
			is_reading_bond = false;
			is_reading_angl = true;
			continue;
		}
		if(each_line.find("!NPHI")!=std::string::npos) {
			max_members = get_psf_section_count(each_line);
			psf_info.nphi = max_members;
			std::cout << "READPSF> Reading PSF DIHEDRAL (" << max_members << ")" << std::endl;
			max_lines = max_members/2 + 1;
			is_reading_angl = false;
			is_reading_dihe = true;
			continue;
		}
		// Now process each line based on flags.
		if(is_reading_atom) parse_PsfAtom(each_line, max_lines);
		if(is_reading_bond) parse_PsfBond(each_line, max_lines);
		if(is_reading_angl) parse_PsfAngl(each_line, max_lines);
		if(is_reading_dihe) parse_PsfDihe(each_line, max_lines);
	}
	inp_file.close();
// debug check
#ifdef _DEBUG
	const int w = 6;	// width of index output
	std::cout << "DEBUG> PSF information" << std::endl;
	for(auto it=atom_topology.cbegin(); it!=atom_topology.cend(); ++it)
		std::cout << std::setw(w) << "Atom:" << std::setw(w) <<it->atom_id << std::setw(w) 
			<< it->atom_type << std::endl;
	for(auto it=bond_topology.cbegin(); it!=bond_topology.cend(); ++it)
		std::cout << std::setw(w) << "Bond: " << it->atom_i << std::setw(w) << it->atom_j 
			<< std::endl;
	for(auto it=angl_topology.cbegin(); it!=angl_topology.cend(); ++it)
		std::cout << std::setw(w) << "Angle: " << it->atom_i << std::setw(w) << it->atom_j 
			<< std::setw(w) << it->atom_k << std::endl;
	for(auto it=dihe_topology.cbegin(); it!=dihe_topology.cend(); ++it)
		std::cout << std::setw(w) << "Dihedral: " << std::setw(w) << it->atom_i << std::setw(w) 
			<< it->atom_j << std::setw(w) << it->atom_k << std::setw(w) 
			<< it->atom_l << std::endl;
#endif
	std::cout << "DEBUG> Size of atom topology array: " << atom_topology.size() << std::endl;
	return true;
}

unsigned int get_psf_section_count(std::string &inp_str) {
	auto npos = inp_str.find_first_not_of("0123456789");
	return std::stol(inp_str.substr(0,npos-1));
}

void parse_PsfAtom(const std::string &inp_data, const unsigned int& max_lines) {
	static unsigned int counter = 0;
	if(counter > max_lines) return;
	PsfAtom psf_atom;
	std::stringstream inp_stream(inp_data);
	inp_stream >> psf_atom.atom_id >> psf_atom.seg_name >> psf_atom.resi_id
		>> psf_atom.resi_name >> psf_atom.atom_name >> psf_atom.atom_type
		>> psf_atom.charge >> psf_atom.mass >> psf_atom.imove 
		>> psf_atom.e_negativity >> psf_atom.hardness;
	atom_topology[psf_atom.atom_id]=psf_atom;
	counter++;
}

void parse_PsfBond(const std::string &inp_data, const unsigned int& max_lines) {
	static unsigned int counter = 0;
	if(counter > max_lines) return;
	std::stringstream inp_stream(inp_data);
	unsigned int n_loop = n_of_words(inp_data)/2;
	n_loop = n_loop==BOND_PER_LINE ? BOND_PER_LINE : n_loop;
	for(unsigned int i=0; i<n_loop; ++i) {
		PsfBond psfbond;
		inp_stream >> psfbond.atom_i >> psfbond.atom_j;
		bond_topology.push_back(psfbond);	
	}
	counter++;
}

void parse_PsfAngl(const std::string &inp_data, const unsigned int& max_lines) {
	static unsigned int counter = 0;
	if(counter > max_lines) return;
	std::stringstream inp_stream(inp_data);
	unsigned int n_loop = n_of_words(inp_data)/3;
	n_loop = n_loop==ANGL_PER_LINE ? ANGL_PER_LINE : n_loop;
	for(unsigned int i=0; i<n_loop; ++i) {
		PsfAngl psfangl;
		inp_stream >> psfangl.atom_i >> psfangl.atom_j >> psfangl.atom_k;
		angl_topology.push_back(psfangl);
	}
	counter++;
}

void parse_PsfDihe(const std::string &inp_data, const unsigned int& max_lines) {
	static unsigned int counter = 0;
	if(counter > max_lines) return;
	std::stringstream inp_stream(inp_data);
	unsigned int n_loop = n_of_words(inp_data)/4;
	n_loop = n_loop==DIHE_PER_LINE ? DIHE_PER_LINE : n_loop;
	for(unsigned int i=0; i<n_loop; ++i) {
		PsfDihe psfdihe;
		inp_stream >> psfdihe.atom_i >> psfdihe.atom_j >> psfdihe.atom_k >> psfdihe.atom_l;
		dihe_topology.push_back(psfdihe);
	}
	counter++;
}

/* interfaces to access PSF data structures. */
PsfInfo get_psf_info() {
	return psf_info;
}

const std::string get_atom_type(const unsigned int query_id) {
	return atom_topology[query_id].atom_type;
}

bool fix_atom (unsigned int index) {
	if(index > atom_topology.size()) {
		std::cout << "ERROR> Requested atom index (" << index 
			<< ") exceeds maxium atom number." << std::endl;
		return false;
	}
	atom_topology[index].imove = 1;
	return true;
}

bool is_movable (unsigned int index) {
	return 0 == atom_topology[index].imove;
}

bool is_good_psf() {
	if(atom_topology.empty()) {
		std::cout << "ERROR> No atom information found in PSF file" << std::endl;
		return false;
	}
	if(bond_topology.empty()) 
		std::cout << "WARNING> No bond information found in PSF file" << std::endl;
	if(angl_topology.empty()) 
		std::cout << "WARNING> No angle information found in PSF file" << std::endl;
	if(dihe_topology.empty()) 
		std::cout << "WARNING> No dihedral information found in PSF file" << std::endl;
	return true;
}

const std::vector<PsfBond>& psf_bond() {
	return bond_topology;
}

const std::vector<PsfAngl>& psf_angl() {
	return angl_topology;
}

const std::vector<PsfDihe>& psf_dihe() {
	return dihe_topology;
}