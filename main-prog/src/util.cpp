#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "../include/util.h"
#include "../include/energy.h"
#include "../include/coor.h"
#include "../include/random.h"

void print_banner() {
    std::cout << "-----------------" << std::endl;
    std::cout << "EZAFM version 1.0" << std::endl;
    std::cout << "-----------------" << std::endl;
}

void debug(std::string msg) {
    std::cout << "DEBUG> " << msg << std::endl;
}

/*check if the program has input file*/
bool has_input(const int inp_num) {
    if(inp_num < 2) {
        std::cout << "ERROR> Missing input file." << std::endl;
        return false;
    } 
    return true;
}

/*check if a string can be converted to integer*/
bool is_integer(std::string inp_str) {
    return inp_str.find_first_not_of("0123456789")==std::string::npos;
}

/* function to determine algebraic sign */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/*check if a string can be converted to float*/
bool is_double(std::string inp_str) {
    char* endptr = nullptr;
    std::strtod(inp_str.c_str(), &endptr);
    if(*endptr!='\0' || endptr==inp_str.c_str()) return false;
    return true;
}

bool is_ignore_case_equal(std::string str1, std::string str2) {
    std::transform(str1.begin(), str1.end(), str1.begin(), ::toupper);
    std::transform(str2.begin(), str2.end(), str2.begin(), ::toupper);
    return str1.compare(str2)==0 ? true : false;
}

bool is_single_word(const std::string &inp_str) {
    int count = 0;
    std::stringstream tmp_stream;
    std::string tmp_holder;
    tmp_stream.str(inp_str);
    while(tmp_stream >> tmp_holder) count++;
    return count==1;
}

unsigned int n_of_words(const std::string &inp_str) {
    int counter = 0;
    std::stringstream tmp_stream(inp_str);
    std::string tmp_holder;
    while(tmp_stream >> tmp_holder) counter++;
    return counter;
}

void skip_n_words(std::stringstream &inp_stream, unsigned int nskip) {
	unsigned int counter = 0;
	std::string dummy;
	while(counter < nskip) {
		inp_stream >> dummy;
		counter++;
	}
}

void print_boundary(unsigned int n) {
    std::string s(n, '*');
    std::cout << s << std::endl;
}
