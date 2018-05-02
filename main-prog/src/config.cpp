#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "../include/config.hpp"
#include "../include/dyna.hpp"
#include "../include/uservar.hpp"
#include "../include/util.hpp"

bool Config::read_config(const std::string& inp_name) {
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open file [" << inp_name << "]" 
                  << std::endl;
        return false;
    } else {
        std::cout << "ReadCFG> Reading configuration file from [" << inp_name << "]" 
                  << std::endl; 
    }
    // File successfully opened, now reading commands.
    std::string each_line; 
    std::string cmd_line; 
    std::string cmd_tmp;
    std::stringstream cmd_stream;
    std::vector<std::string> cmds;
    // Skip empty lines, commented lines; Concatenate "-" connected lines.
    while(std::getline(inp_file, each_line)) {
        //remove off-command-line comments.
        if(each_line[0]=='!' || each_line.empty()) continue; 
        size_t newlen = each_line.find_first_of('!');
        // remove in-command-line comments.
        cmd_line = each_line.substr(0,newlen);
        cmd_stream.clear();
        cmd_stream.str(cmd_line);
        while(cmd_stream >> cmd_tmp)
            cmds.push_back(cmd_tmp);   
        if(cmds.back()=="-")
            continue;
        // remove "-" if applicable.
        cmds.erase(std::remove(cmds.begin(),cmds.end(),"-"),cmds.end()); 
		_config.push_back(cmds);
        cmds.clear();
    }
    inp_file.close();
    return true;
}

bool Config::exec_config(void) {
    DynaSystem dyna_system;
    UserVar  user_var;
    /* Now parse commands set*/
    for (auto& each_cmd : this->_config) {
        auto cmd = each_cmd.front();
        if( "set"==cmd ) {
            user_var.update(each_cmd);
        } else if ("system"==cmd) {
            std::string inp_type=each_cmd[1]; 
            std::string inp_name=each_cmd[2]; 
            inp_name=user_var.query(inp_name);
            if("psf"==inp_type) dyna_system.read_psf(inp_name);
            if("cor"==inp_type) dyna_system.read_cor(inp_name);
            if("prm"==inp_type) dyna_system.read_prm(inp_name);
            if("fix"==inp_type) {
                size_t fix_id = std::stol(inp_name);
                dyna_system.set_fix_atom(fix_id);
            }
        } else if ("afm"==cmd) {
            dyna_system.setup_afm(each_cmd, user_var);
        } else if ("dyna"==cmd) {
            dyna_system.setup_dyna(each_cmd, user_var);
            dyna_system.run_dyna();
        } else if ("stop"==cmd) {
            return true;
        } else {
            std::cout << "ERROR> Unknown command (" << cmd << ")"
                       << std::endl;
            return false;
        }
    }
    return true;
}