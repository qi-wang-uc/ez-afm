#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "../include/read.h"
#include "../include/util.h"
#include "../include/uservar.h"
#include "../include/param.h"
#include "../include/psf.h"
#include "../include/coor.h"

bool read_data(std::vector<std::string> cmds) {
    std::string inp_form = cmds[1];
    std::string inp_name = get_uservar(cmds[2]);
    if (inp_form=="prm") {
        if(!read_prm(inp_name)) return false;
    } else if (inp_form=="psf") {
        if(!read_psf(inp_name)) return false;
    } else if (inp_form=="cor") {
        if(!read_cor(inp_name)) return false;
    } else {
        std::cout << "Unrecorgnized or unnecessary file type!" << std::endl;
        return false;
    }
    return true;
}

/*********************************************************************
   Note: since the PSF file contains all the information of RTF file, 
   it seems un-necessary to read RTF in the program. However, one can 
   implement the function here for other purposes.
*********************************************************************/