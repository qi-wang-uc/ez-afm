#ifndef DCD_H
#define DCD_H

#include <string>
#include <fstream>

/* dcd_trajectory_file = dcd_header + dcd_pads + dcd_coordinates */

struct DCD_Header {
    unsigned int natom;
    unsigned int nfile;
    unsigned int npriv;
    unsigned int nsavc;
    unsigned int nstep;
    int ifpbc;
    double delta;    // delta =  dt/TIMEFACT
    std::string prog_title; // need to resize to 80.
    std::string user_title; // need to resize to 80.
    /*
    DCD_Header() {}

    DCD_Header(int natom, int nfile, int npriv, int nsavc, int nstep, int ifpbc, 
        double delta, std::string prog_title, std::string user_title) :
        natom(natom),nfile(nfile),npriv(npriv),nsavc(nsavc),nstep(nstep),ifpbc(ifpbc),
        delta(delta),prog_title(prog_title),user_title(user_title)
        {}
    */
};

struct DCD_Pads {
    const int pad0   = 0;
    const int pad2   = 2;
    const int pad4   = 4;
    const int pad24  = 24;
    const int pad84  = 84;
    const int pad164 = 164;  
};

void write_dcdheader(std::ofstream& dcd_file, DCD_Header& dcd_header, DCD_Pads& dcd_pads);

void write_dcdframe(std::ofstream& dcd_file, unsigned int natom, DCD_Pads& dcd_pads);

#endif