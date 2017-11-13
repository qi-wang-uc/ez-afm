#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "../include/random.h"

/*******************************************************************
               GAUSSIAN RANDOM NUMBER DATA 
*******************************************************************/
static std::vector<real> gaussian_rand;

void init_random(const unsigned int natom) {
    std::cout << "RANDGEN> Initialising random number array." << std::endl;
    if(!gaussian_rand.empty()) return;
    gaussian_rand.resize(3*natom);
    std::cout << "RANDGEN> Done." << std::endl << std::endl;
}

void gen_random(real mean, real dev) {
    if(gaussian_rand.empty()) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<real> nd(mean, dev);
    for(unsigned int i=0; i<gaussian_rand.size(); ++i) {
        gaussian_rand[i]=nd(gen);
    }
}

const real& get_random(unsigned int query) {
    return gaussian_rand[query];
}

void print_random() {
    std::cout << "*** Printing random number array ***" << std::endl;
    for(auto it=gaussian_rand.cbegin(); it!=gaussian_rand.cend(); ++it)
        std::cout << *it << std::endl;
    std::cout << "************************************" << std::endl;

}