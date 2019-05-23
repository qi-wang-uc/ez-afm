#ifndef USERVAR_HPP
#define USERVAR_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "define.hpp"

class UserVar {
    private:
        StrTable _table;
    public:
        void print(void);
        void update(StrVec cmds);
        Str  query(Str& name) const;      
};

#endif
