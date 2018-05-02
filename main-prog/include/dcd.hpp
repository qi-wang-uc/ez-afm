#ifndef DCD_HPP
#define DCD_HPP

#include <string>
#include <fstream>

/* dcd_trajectory_file = dcd_header + dcd_pads + dcd_coordinates */

struct DCD_Header {
    int32_t natom;
    int32_t nfile;
    int32_t npriv;
    int32_t nsavc;
    int32_t nstep;
    int32_t ifpbc;
    float delta; // delta =  dt/TIMEFACT
    char* prog_title; 
    char* user_title; 
    DCD_Header() {}
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
    protected:
        DCD_Header _dcd_header;
        DCD_Pads   _dcd_pads;
    public:
        void write_dcdheader(std::ofstream& dcd_file); 
        void write_dcdframe(std::ofstream& dcd_file, const float* fX, const float* fY, const float* fZ, const size_t& natom);
};

#endif