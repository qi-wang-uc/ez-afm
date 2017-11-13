#include <iostream>
#include <fstream>
#include <string>
#include "../include/dcd.h"
#include "../include/coor.h"

void write_dcdheader(std::ofstream& dcd_file, DCD_Header& dcd_header, DCD_Pads& dcd_pads) {
    /* dcd header part 1. */
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad84), 4);
    dcd_file.write("CORD", 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.nfile), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.npriv), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.nsavc), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.nstep), 4);
    for(int i=0; i<5; ++i)
        dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad0), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.delta), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.ifpbc), 4);
    for(int i=0; i<8; ++i)
        dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad0), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad24), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad84), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad164), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad2), 4);
    /* dcd header part 2. */
    dcd_file.write(dcd_header.prog_title.c_str(), 80);
    dcd_file.write(dcd_header.user_title.c_str(), 80);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad164), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_header.natom), 4);
    dcd_file.write(reinterpret_cast<const char*>(&dcd_pads.pad4), 4);
}

void write_dcdframe(std::ofstream& dcd_file, unsigned int natom, DCD_Pads& dcd_pads) {
    const int pad4N = 4*natom;
    /* write x coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(p_xcoor()), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    /* write y coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(p_ycoor()), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    /* write z coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(p_zcoor()), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
}
