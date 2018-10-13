#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>
#include "define.hpp"

using CMDs = std::vector<std::vector<std::string> >;

class Config {
    private:    
        CMDs _config;
    public:
        bool read_config(const Str& inp_name);
        bool exec_config(void);
};

#endif
