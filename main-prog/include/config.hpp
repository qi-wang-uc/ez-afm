#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

using CMDs = std::vector<std::vector<std::string> >;

class Config {
    private:    
        CMDs _config;
    public:
        bool read_config(const std::string& inp_name);
        bool exec_config(void);
};

#endif
