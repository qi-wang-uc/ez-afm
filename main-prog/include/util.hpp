#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <ctime>
#include "define.hpp"

using TimeStamp = std::chrono::system_clock;
using SysClock  = std::chrono::time_point<std::chrono::system_clock>;

bool is_integer(Str inp_str);

bool is_double(Str inp_str);

bool is_ignore_case_equal(Str str1, Str str2);

bool is_single_word(const Str &inp_str);

Int n_of_words(const Str &inp_str);

void skip_n_words(std::stringstream &inp_stream, Int nskip);

void print_boundary(Int n);

template <typename T> 
Int sgn(T val);

void time_elapsed(SysClock tick,  SysClock tock);

Str time_stamp();

void debug(Int code);

Str int2str(Int integer);

#endif