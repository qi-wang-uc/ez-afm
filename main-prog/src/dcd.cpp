#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../include/dcd.hpp"

void DcdData::write_dcdheader(std::ofstream& dcd_file) {
    /* dcd header part 1. */
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad84), 4);
    dcd_file.write("CORD", 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.nfile), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.npriv), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.nsavc), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.nstep), 4);
    for(int i=0; i<5; ++i)
        dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad0), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.delta), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.ifpbc), 4);
    for(int i=0; i<8; ++i)
        dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad0), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad24), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad84), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad164), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad2), 4);
    /* dcd header part 2. */
    dcd_file.write(this->_dcd_header.prog_title, 80);
    dcd_file.write(this->_dcd_header.user_title, 80);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad164), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad4), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_header.natom), 4);
    dcd_file.write(reinterpret_cast<const char*>(&this->_dcd_pads.pad4), 4);
}

void DcdData::write_dcdframe(std::ofstream& dcd_file, const float* fX, const float* fY, const float* fZ, const size_t& natom) {
    const int32_t pad4N = 4*natom;
    /* write x coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(fX), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    /* write y coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(fY), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    /* write z coordinates. */
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
    dcd_file.write(reinterpret_cast<const char*>(fZ), pad4N);
    dcd_file.write(reinterpret_cast<const char*>(&pad4N), 4);
}
