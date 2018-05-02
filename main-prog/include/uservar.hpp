#ifndef USERVAR_HPP
#define USERVAR_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>

// Not supposed to be inherited.
class UserVar {
    private:
        std::map<std::string,std::string> _table;
    public:
        void print(void);
        void update(std::vector<std::string> cmds);
        std::string query(std::string& name) const;
        
};

#endif
