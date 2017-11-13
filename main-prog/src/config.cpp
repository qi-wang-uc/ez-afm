#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "../include/config.h"

std::vector<std::vector<std::string> > read_config(const std::string inp_name) {
	std::vector<std::vector<std::string> > result;
    std::ifstream inp_file(inp_name);
    if(!inp_file.is_open()) {
        std::cout << "ERROR> Cannot open file [" << inp_name << "]" << std::endl;
        return result;
    } else {
        std::cout << "READCFG> Reading configuration file from [" << inp_name << "]" << std::endl; 
    }
    // File successfully opened, now reading commands.
    std::string each_line; 
    std::string cmd_line; 
    std::string cmd_tmp;
    std::stringstream cmd_stream;
    std::vector<std::string> cmds;
    // Skip empty lines, commented lines; Concatenate "-" connected lines.
    while(std::getline(inp_file, each_line)) {
        if(each_line[0]=='!' || each_line.empty()) continue; //remove off-command-line comments.
        int newlen = each_line.find_first_of('!');
        cmd_line = each_line.substr(0,newlen);  // remove in-command-line comments.
        cmd_stream.clear();
        cmd_stream.str(cmd_line);
        while(cmd_stream >> cmd_tmp) cmds.push_back(cmd_tmp);   // store cmd in vector
        if(cmds.back()=="-") continue;  // if continued line continue reading next line.
        cmds.erase(std::remove(cmds.begin(),cmds.end(),"-"),cmds.end()); // remove "-" if has.
		result.push_back(cmds);
        cmds.clear();	// clear for next command parsing
    }
    inp_file.close();
	return result;
}
