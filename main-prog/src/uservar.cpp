#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include "../include/uservar.hpp"

void UserVar::print(void) {
    std::cout << this->_table.size() << std::endl;
	std::cout << "EZAFM> User-defined variables " << std::endl;
	for (auto& uv : this->_table) {
        std::cout << std::setw(10) << uv.first 
                  << " : " 
                  << std::setw(10) << uv.second
                  << std::endl;
    }
    std::cout << "EZAFM> " << std::endl;
}

void UserVar::update(std::vector<std::string> cmds) {
	this->_table[cmds[1]] = ("="==cmds[2]) ? cmds[3]:cmds[2];
}

std::string UserVar::query(std::string& name) const {
    std::string invalid_query = "";
    std::cout << "UserVar> Retriving (" << name << ")" << std::endl;
    while(name.find("@")!=std::string::npos) {
        auto var_begin = name.find_first_of("@");
        auto val_begin = var_begin + 2;
        auto left_brace  = name.find_first_of("{");
        auto right_brace = name.find_first_of("}");
        auto len_of_val = right_brace - left_brace - 1;
        auto len_of_var = len_of_val + 3;
        std::string var = name.substr(var_begin, len_of_var);
        std::string val = name.substr(val_begin, len_of_val);
        if(this->_table.find(val)==this->_table.end()) {
            std::cout << "ERROR> Cannot find (" << val 
                      << ") in user defined variables." 
                      << std::endl;
            return invalid_query;
        }
        name = name.replace(name.find_first_of("@"), len_of_var, this->_table.at(val));
    }
    return name;
}

