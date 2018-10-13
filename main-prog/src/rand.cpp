#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>
#include "../include/rand.hpp"

void Rand::init_rand(const Int& natom) {
    std::cout << "RandGen> Initialising random number array." << std::endl;
    if(!(_xrand.empty()&&_yrand.empty()&&_zrand.empty())) return;
    _xrand.resize(natom);
    _yrand.resize(natom);
    _zrand.resize(natom);
    _size = natom;
    std::cout << "RandGen> Done." << std::endl;
}

void Rand::gen_rand(const Real& mean, const Real& dev) {
    if(_xrand.empty()||_yrand.empty()||_zrand.empty()) return;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<Real> nd(mean, dev);
    std::generate(_xrand.begin(), _xrand.end(), [&](){return nd(gen);});
    std::generate(_yrand.begin(), _yrand.end(), [&](){return nd(gen);});
    std::generate(_zrand.begin(), _zrand.end(), [&](){return nd(gen);});
}

Vec3d Rand::get_rand(const Int& atomid) const {
    return Vec3d(_xrand.at(atomid), _yrand.at(atomid), _zrand.at(atomid));
}

void Rand::print_rand(void) const {
    std::cout << "RandGen> Printing random number array " << std::endl;
    for(Int i = 0; i<_size; ++i)
        std::cout << i << std::endl;
    std::cout << "************************************" << std::endl;

}

const Real* Rand::px() const {
    return this->_xrand.data();
};
const Real* Rand::py() const {
    return this->_yrand.data();
};
const Real* Rand::pz() const {
    return this->_zrand.data();
};