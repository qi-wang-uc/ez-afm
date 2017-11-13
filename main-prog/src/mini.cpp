#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include "../include/mini.h"
#include "../include/sd.h"

void prep_mini(std::vector<std::string> cmds) {
    std::string method = cmds[1];
    unsigned int steps = std::stol(cmds[2]);
    run_mini(method, steps);
}

void run_mini(std::string method, const unsigned int nstep) {
    if(method=="sd") {
        std::cout << "EMINI> Energy minimization using [Steepest Descent] method for " 
         << nstep << " steps." << std::endl;
         steepest_descent(nstep);
    } else if (method=="abnr") {
        std::cout << "EMINI> Energy minimization using [Adopted Basis Newton-Raphson ] method for " 
         << nstep << " steps." << std::endl;
         //newton_raphson();
    } else if (method=="conj") {
        std::cout << "EMINI> Energy minimization using [Conjugate Gradient] method for " 
         << nstep << " steps." << std::endl;
         //conjugate_gradient();
    } else {
        std::cout << "ERROR> Unrecorgnized minimization method." << std::endl;
    }
}
