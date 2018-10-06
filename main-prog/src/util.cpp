#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include "../include/util.hpp"

/*check if a string can be converted to integer*/
bool is_integer(Str inp_str) {
    return inp_str.find_first_not_of("0123456789")==std::string::npos;
}

/* function to determine algebraic sign */
template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/*check if a string can be converted to float*/
bool is_double(Str inp_str) {
    char* endptr = nullptr;
    std::strtod(inp_str.c_str(), &endptr);
    if(*endptr!='\0' || endptr==inp_str.c_str()) return false;
    return true;
}

bool is_ignore_case_equal(Str str1, Str str2) {
    std::transform(str1.begin(), str1.end(), str1.begin(), ::toupper);
    std::transform(str2.begin(), str2.end(), str2.begin(), ::toupper);
    return str1.compare(str2)==0 ? true : false;
}

bool is_single_word(const Str &inp_str) {
    int count = 0;
    std::stringstream tmp_stream;
    Str tmp_holder;
    tmp_stream.str(inp_str);
    while(tmp_stream >> tmp_holder) count++;
    return count==1;
}

Int n_of_words(const Str &inp_str) {
    int counter = 0;
    std::stringstream tmp_stream(inp_str);
    Str tmp_holder;
    while(tmp_stream >> tmp_holder) counter++;
    return counter;
}

void skip_n_words(std::stringstream &inp_stream, Int nskip) {
	Int counter = 0;
	Str dummy;
	while(counter < nskip) {
		inp_stream >> dummy;
		counter++;
	}
}

void print_boundary(Int n) {
    Str s(n, '*');
    std::cout << s << std::endl;
}

void time_elapsed(SysClock tick, SysClock tock) {
    auto duration = tock - tick;
    auto seconds =  std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto total_seconds = seconds.count();
    auto n_hours   =  (total_seconds / 3600);
    auto n_minutes =  (total_seconds - n_hours*3600)/60;
    auto n_seconds =  (total_seconds - n_hours*3600 - n_minutes*60);
    std::cout << "EZAFM> Time elapsed: "
              << "(" << n_hours   << ") HR : "
              << "(" << n_minutes << ") MIN : "
              << "(" << n_seconds << ") SEC" << std::endl;
}

void debug(Int code) {
    std::cout << ">>>>>>>> DEBUG :" << code << std::endl;
}

Str time_stamp() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time), "%Y-%m-%d %X");
    return ss.str();
}

Str int2str(Int integer) {
    return "(" + std::to_string(integer) + ")";
}
