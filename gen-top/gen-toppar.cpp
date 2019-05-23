#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include "gen-toppar.hpp"

/*********************************************************************
	Note: the format are derived from exist CHARMM-compatible 
	rtf/prm/psf/crd files thus cannot be used as reference.
*********************************************************************/

/****************** USER SPECIFIED VARABILES ************************/
std::string pdb_name = "1ubq.pdb";	// make sure to use .pdb as extension.
const double eps_native = -2.0;		// well-depth of native potential. must be NEGATIVE.
const double r_cutoff = 8.0;		// cutoff value to determine native contacts.
/********************************************************************/

std::vector<amino> amino_acids;  // stores all the CAs in pdb file.

int main(int argc, char* argv[]) {
    if(!read_pdb(pdb_name)) return -1; 
    std::string job_name = get_jobname(pdb_name);

    make_cor(job_name);
    make_psf(job_name);
    make_rtf(job_name);
    make_prm(job_name);

    return 0;
}

double calc_dist(const amino& a1, const amino& a2) {
    double dx = a2.x - a1.x;
    double dy = a2.y - a1.y;
    double dz = a2.z - a1.z;
    double r2 = dx*dx + dy*dy + dz*dz;
    return sqrt(r2);
}

bool read_pdb(const std::string& pdb_name) {
    std::ifstream pdb_file(pdb_name);
    if(!pdb_file) {
        std::cerr << "Could not open pdf file [" << pdb_name << "]" << std::endl;
        return false;
    }
    std::string each_line;
    unsigned int i = 1;
    std::string s = "PROT";
    double w = 0.00000;
    while(getline(pdb_file, each_line)) {
        if(each_line.substr(0,4)=="ATOM" && each_line.substr(13,2)=="CA") {
            float x = atof(each_line.substr(31,8).c_str());
            float y = atof(each_line.substr(39,8).c_str());
            float z = atof(each_line.substr(47,8).c_str());
            std::string res = each_line.substr(17,3).c_str();
            int resid = strtol(each_line.substr(23,3).c_str(),NULL,10);
            amino_acids.push_back(amino(i,i,res,"CA",x,y,z,s,resid,w));
            i++;
        }
    }
    pdb_file.close();
    auto natom = i-1;
    std::cout << "After reading pdb file, (" << natom << ") atoms were found." << std::endl;
	return true;
}

std::string get_jobname(const std::string& pdb_name) {
    size_t npos = pdb_name.find(".pdb");
    return (npos==std::string::npos) ? pdb_name : pdb_name.substr(0,npos);
}

std::string get_header(std::string job_name) {
    time_t now = time(0);
    std::string dt = ctime(&now);
    dt.pop_back();
    std::transform(job_name.begin(), job_name.end(), job_name.begin(), ::toupper);
    std::string head_info = std::string("* FILES FOR JOB: ") + job_name;
    std::string time_info = std::string("* DATE: ") + dt + fmt_space.space4;
    std::string user_info = "CREATED BY USER: " + std::string(getenv("USER"));
    return head_info + std::string("\n") + time_info + user_info + std::string("\n");
}

void make_cor(const std::string& job_name) {
    std::string out_name = job_name + ".cor";
    std::ofstream out_file(out_name);
    out_file << get_header(job_name) << "*" << std::endl;
    out_file << std::setw(5) << amino_acids.size() << std::endl;
    out_file << std::setiosflags(std::ios::fixed);
    for(int i=0; i<amino_acids.size(); i++){
        out_file << std::right << std::setw(5) << amino_acids[i].atomno
                 << std::right << std::setw(5) << amino_acids[i].resno
                 << fmt_space.space1 << std::left << std::setw(4) << amino_acids[i].res
                 << fmt_space.space1 << std::left << std::setw(4) << amino_acids[i].type
                 << std::right << std::setw(10) << std::setprecision(5) << amino_acids[i].x
                 << std::right << std::setw(10) << std::setprecision(5) << amino_acids[i].y
                 << std::right << std::setw(10) << std::setprecision(5) << amino_acids[i].z
                 << fmt_space.space1 << std::setw(4) << amino_acids[i].segid
                 << fmt_space.space1 << std::left << std::setw(4) << amino_acids[i].resid
                 << std::right << std::setw(10) << std::setprecision(5) << amino_acids[i].w
                 << std::endl;
    }
    out_file.close();
}

void make_psf(const std::string& job_name) {
    std::string out_name = job_name + ".psf";
    std::ofstream out_file(out_name);
//HEADER
    out_file << "PSF CMAP CHEQ XPLOR" << std::endl << std::endl;
    out_file << std::setw(8) << 2 << " !NTITLE" << std::endl;
    out_file << get_header(job_name) << std::endl;
//NATOM
    out_file << std::setw(8) << amino_acids.size() << " !NATOM" << std::endl;
// CHARMM PSF VARIABLE 
    int imove_i = 0;  // imove=0 : free atom
    int atc_i   = 0;  // atc: atom type code
    double cg_i = 0.00000;
    double amass_i = 100.00;
    double ech_i = 0.00000;  // electronegativity of atoms.
    double eha_i = -0.301140E-02;  // hardness for atoms.
    for(int i=0; i<amino_acids.size(); i++) {
        atc_i = i+1;
        out_file << std::right << std::setw(8) << amino_acids[i].atomno
                 << fmt_space.space1 << std::setw(4) << amino_acids[i].segid
                 << fmt_space.space1 << std::left << std::setw(4) << amino_acids[i].resid
                 << fmt_space.space1 << std::left << std::setw(4) << amino_acids[i].res
                 << fmt_space.space1 << std::left << std::setw(4) << amino_acids[i].type
                 << fmt_space.space1 << std::left << std::setw(4) << std::string("A")+std::to_string(atc_i)
                 << "    0.00000       100.000           0   0.00000     -0.301140E-02" // quick & dirty
                 << std::endl;
    }
    out_file << std::endl;
	out_file << std::right;
//NBOND
    out_file << std::setw(8) << amino_acids.size()-1 << " !NBOND: bonds" << std::endl;
    for(unsigned int i=0; i<amino_acids.size()-1; i++) {
        unsigned int ibond = i + 1;
        unsigned int jbond = i + 2;
        out_file << std::setw(8) << ibond << std::setw(8) << jbond;
        if(ibond%4==0) out_file << std::endl;
    }
    out_file << std::endl << std::endl;
//NTHETA
    out_file << std::setw(8) << amino_acids.size()-2 << " !NTHETA: angles" << std::endl;
    for(unsigned int i=0; i<amino_acids.size()-2; i++) {
        unsigned int itheta = i+1;
        unsigned int jtheta = i+2;
        unsigned int ktheta = i+3;
        out_file << std::setw(8) << itheta << std::setw(8) << jtheta << std::setw(8) << ktheta;
        if(itheta%3==0) out_file << std::endl;
    }
    out_file << std::endl << std::endl;
//NPHI
    out_file << std::setw(8) << amino_acids.size()-3 << " !NPHI: dihedrals" << std::endl;
    for(unsigned int i=0; i<amino_acids.size()-3; i++) {
        unsigned int iphi = i+1;
        unsigned int jphi = i+2;
        unsigned int kphi = i+3;
        unsigned int lphi = i+4;
        out_file << std::setw(8) << iphi << std::setw(8) << jphi << std::setw(8) << kphi << std::setw(8) << lphi;
        if(iphi%2==0) out_file << std::endl;
    }
    out_file << std::endl << std::endl;
//NIMPHI
    out_file << std::setw(8) << 0 << " !NIMPHI: impropers" << std::endl;
    out_file << std::endl << std::endl;
//NDON
    out_file << std::setw(8) << 0 << " !NDON: donors" << std::endl;
    out_file << std::endl << std::endl;
//NACC
    out_file << std::setw(8) << 0 << " !NACC: acceptors" << std::endl;
    out_file << std::endl << std::endl;
//NNB
    out_file << std::setw(8) << 0 << " !NNB" << std::endl;
    out_file << std::endl;
//IBLO (Nonbonded exclusion)
    for(unsigned int i=0; i<amino_acids.size(); i++) {
        out_file << std::setw(8) << 0;
        if((i+1)%8==0) out_file << std::endl;
    }
//NGRP NST2
    out_file << std::endl << std::endl;
    out_file << std::setw(8) << amino_acids.size() << std::setw(8) << 0 << " !NGRP NST2" << std::endl;
    for(unsigned int i=0; i< amino_acids.size(); i++) {
        out_file << std::setw(8) << i << std::setw(8) << 0 << std::setw(8) << 0;
        if((i+1)%3==0) out_file << std::endl;
    }
//MOLNT
    out_file << std::endl << std::endl;
    out_file << std::setw(8) << 1 << " !MOLNT" << std::endl;
    for(unsigned int i=0; i<amino_acids.size(); i++) {
        out_file << std::setw(8) << 1;
        if((i+1)%8==0) out_file << std::endl;
    }
    out_file << std::endl << std::endl;
//NUMLP NUMLPH
    out_file << std::setw(8) << 0 << std::setw(8) << 0 << " !NUMLP NUMLPH" << std::endl;
    out_file << std::endl;
    out_file << std::setw(8) << 0 << " !NCRTERM: cross-terms" << std::endl;
    out_file << std::endl;
    out_file.close();
}

void make_rtf(const std::string& job_name) {
    std::string out_name = job_name + ".rtf";
    std::ofstream out_file(out_name);
    std::string title_name = job_name;
    std::transform(title_name.begin(), title_name.end(), title_name.begin(), ::toupper);
    out_file << "* CHARMM TOPOLOGY FILE FOR " << title_name << std::endl << "*" << std::endl;
    out_file << "   20   1" << std::endl;
    for(int i=0; i<amino_acids.size(); i++) {
        out_file << "MASS " << std::right << std::setw(4) << i+1 << fmt_space.space1
                 << std::left << std::setw(4) << std::string( "A") + std::to_string(i+1) 
                 << fmt_space.space5 << "100.0"
                 << std::endl;
    }
    out_file << std::endl << "DECL +CA" << std::endl;
    out_file << std::endl << "AUTOGENERATE ANGLES DIHEDRAL" << std::endl << std::endl;

    for(int i=0; i<amino_acids.size(); i++){
        out_file << "RESI " << std::string("P") + std::to_string(i+1) << fmt_space.space1 
                 << "0.00" << std::endl;
        out_file << "GROU" << std::endl;
        out_file << "ATOM " << "CA" << fmt_space.space1 << std::string("A") + std::to_string(i+1) 
                 << " 0.00" << std::endl;
        out_file << "BOND " << "CA +CA" << std::endl;
        out_file << std::endl;
    }
    out_file.close();
}

void make_prm(const std::string& job_name) {
    std::string out_name = job_name + ".prm";
    std::ofstream out_file(out_name);
    std::string title_name = job_name;
    std::transform(title_name.begin(), title_name.end(), title_name.begin(), ::toupper);
    out_file << "* CHARMM PARAM FILE FOR " << title_name << std::endl << "*" << std::endl << std::endl;
//BOND
    out_file << "BOND" << std::endl;
    for(int i=0; i<amino_acids.size()-1; i++) {
        out_file << std::left; 
        out_file << std::setw(4) << std::string("A") + std::to_string(i+1) << fmt_space.space1 
                 << std::setw(4) << std::string("A") + std::to_string(i+2) << "   100.0    3.8" 
                 << std::endl;
    }
    out_file << std::endl;
//ANGLE
    out_file << "ANGLE" << std::endl;
    for(unsigned int i=0; i<amino_acids.size()-2; i++) {
        out_file << std::left;
        out_file << std::setw(4) << std::string("A") + std::to_string(i+1) << fmt_space.space1 
                 << std::setw(4) << std::string("A") + std::to_string(i+2) << fmt_space.space1 
                 << std::setw(4) << std::string("A") + std::to_string(i+3) << fmt_space.space1 
                 << "   12.5     1.8326"	// 1.8326 in Rad or 105 in Degrees.
                 << std::endl;
    }
    out_file << std::endl;
//DIHEDRAL
    out_file << "DIHEDRAL" << std::endl;
    for(unsigned int i=0; i<amino_acids.size()-3; i++) {
        out_file << std::setw(4) << std::string("A") + std::to_string(i+1) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+2) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+3) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+4) + fmt_space.space1
                 << "  0.90 3  -1.16"  << std::endl ;   //-1.16 in Rad or -66.5 in Degrees
        out_file << std::setw(4) << std::string("A") + std::to_string(i+1) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+2) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+3) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+4) + fmt_space.space1
                 << "  2.27 2  -1.19"  << std::endl ;   //-1.19 in Rad or -68.1 in Degrees
        out_file << std::setw(4) << std::string("A") + std::to_string(i+1) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+2) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+3) + fmt_space.space1
                 << std::setw(4) << std::string("A") + std::to_string(i+4) + fmt_space.space1
                 << "  2.91 1  -0.62"  << std::endl ;   //-0.62 in Rad or -37.3 in Degrees
    }
    out_file << std::endl;
//NONBONDED
    out_file << "NBONDED" << std::endl;
    for(int i=0; i<amino_acids.size(); i++) {
        out_file << std::left;
        out_file << std::setw(4) << std::string("A") + std::to_string(i+1) 
                 << "   0.0    -1e-12     22.71" << std::endl;
    }
    out_file << std::endl;
//NBFIX
    out_file << "NBFIX" << std::endl;
    for(unsigned int i=0; i<amino_acids.size()-3; i++) {
        for(unsigned int j=i+3; j<amino_acids.size(); j++) {
            auto dist = calc_dist(amino_acids[i], amino_acids[j]);
            if(dist <= r_cutoff) {
            out_file << std::setw(4) << std::string("A") + std::to_string(i+1) << fmt_space.space1 
                     << std::setw(4) << std::string("A") + std::to_string(j+1) << fmt_space.space2 
                     << std::to_string(eps_native) << fmt_space.space2 << std::to_string(dist) 
                     << std::endl;
                }
            }
        }
    out_file << std::endl;
//DONE
    out_file << "END" << std::endl;
    out_file.close();
}
