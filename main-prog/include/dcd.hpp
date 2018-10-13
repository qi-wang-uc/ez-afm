#ifndef DCD_HPP
#define DCD_HPP

#include <string>
#include <fstream>
#include "define.hpp"

/* dcd_trajectory_file = dcd_header + dcd_pads + dcd_coordinates */

struct DCD_Header {
    int32_t natom;
    int32_t nfile;
    int32_t npriv;
    int32_t nsavc;
    int32_t nstep;
    int32_t ifpbc;
    float delta; // delta =  dt/TIMEFACT
    const char* prog_title; 
    const char* user_title; 
    DCD_Header() {}
    DCD_Header(int32_t natom, int32_t nfile, int32_t npriv, int32_t nsavc, int32_t nstep, int32_t ifpbc, 
                    float delta, const char* prog_title, const char* user_title ):
                natom(natom),nfile(nfile),npriv(npriv),nsavc(nsavc),nstep(nstep),ifpbc(ifpbc),
                    delta(delta),prog_title(prog_title), user_title(user_title) {}
};

struct DCD_Pads {
    const int32_t pad0   = 0;
    const int32_t pad2   = 2;
    const int32_t pad4   = 4;
    const int32_t pad24  = 24;
    const int32_t pad84  = 84;
    const int32_t pad164 = 164;  
};

class DcdData {
    private:
        DCD_Header _dcd_header;
        DCD_Pads   _dcd_pads;
    public:
        void set_dcdheader(DCD_Header&& inp_header);
        void write_dcdheader(std::ofstream& dcd_file); 
        void write_dcdframe(std::ofstream& dcd_file, const float* fX, const float* fY, const float* fZ, const Int& natom);
};

#endif