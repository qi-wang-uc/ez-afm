#ifndef GENTOPPAR_HPP
#define GENTOPPAR_HPP

#include <string>

/***************** Under development ********************************/
/* string pdb_type = "rcsb";    
   Specify pdb types to read chain names. Currently the program can
   still read both types but only residue type and coordinates are
   extracted from pdb file. For "charmm" pdb types, the bln_table 
   should specifiy 3 histodines. */
/********************************************************************/

struct amino {
    int atomno;         // 1, 2, ..., n.
    int resno;          // 1, 2, ..., n.
    std::string res;    // A1, A2, ..., An.
    std::string type;   // type is CA for all beads.
    double x;           // Cartesian coordinates
    double y;           // ...
    double z;           // ...
    std::string segid;  // use [PROT] by default.
    int resid;          // 1, 2, ..., n.
    double w;           // 0.00000 by default.
    amino(int atomno, int resno, std::string res, std::string type,
        double x, double y, double z, std::string segid, int resid, double w):
        atomno(atomno), resno(resno), res(res), type(type), x(x), y(y), z(z),
        segid(segid), resid(resid), w(w) {}
};

const struct Format_Space {
    std::string space1 = " ";
    std::string space2 = "  ";
    std::string space3 = "   ";
    std::string space4 = "    ";
    std::string space5 = "     ";
    std::string space6 = "      ";
} fmt_space;

double calc_dist(const amino& a1, const amino& a2);

bool read_pdb(const std::string& pdb_name);
std::string  get_jobname(const std::string& pdb_name);
std::string  get_header(std::string job_name);

void make_cor(const std::string& job_name);
void make_psf(const std::string& job_name);
void make_rtf(const std::string& job_name);
void make_prm(const std::string& job_name);

#endif
