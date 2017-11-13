#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include "../include/uservar.h"

/* User defined variable table */ 
static std::map<std::string,std::string> uservar;

/* Add user defined variable */
void set_uservar(std::vector<std::string> cmds) {
	uservar[cmds[1]] = ("="==cmds[2]) ? cmds[3]:cmds[2];
}

/* Variable substitution */
std::string get_uservar(std::string name) {
    std::string false_entry = "";
    while(name.find("@")!=std::string::npos) {
        /*  The usage of user-defined variables should be like:
            xxx@{var1}xxx@{var2}...
            In every iteration, the firt occurence of '@' would be an identifier as 
            beginning and '}' as the end of a variable.
            Here we use [ var = @{val} ] for variable nomenclature.
        */
        auto var_begin = name.find_first_of("@");
        auto val_begin = var_begin + 2;
        auto left_bracket  = name.find_first_of("{");
        auto right_bracket = name.find_first_of("}");
        auto len_of_val = right_bracket - left_bracket - 1;
        auto len_of_var = len_of_val + 3;
        std::string var = name.substr(var_begin, len_of_var);
        std::string val = name.substr(val_begin, len_of_val);
        if(uservar.find(val)==uservar.end()) {
            std::cout << "ERROR> Cannot find (" << val << ") in user defined variables." << std::endl; 
            return false_entry;
        }
        name = name.replace(name.find_first_of("@"), len_of_var, uservar[val]);
    }
    return name;
}

void print_uservar() {
	std::cout << "*** User-defined variables ***" << std::endl;
	for (auto it_uv=uservar.begin(); it_uv!=uservar.end(); ++it_uv)
		std::cout << std::setw(10) << it_uv->first << " : " 
			<<  it_uv->second << std::setw(10) << std::endl;
    std::cout << "******************************" << std::endl;
}
