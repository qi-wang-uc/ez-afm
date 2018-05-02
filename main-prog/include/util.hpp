#ifndef UTIL_HPP
#define UTIL_HPP

#include <vector>
#include <string>
#include <map>
#include <chrono>

using TimeStamp = std::chrono::system_clock;
using SysClock  = std::chrono::time_point<std::chrono::system_clock>;

bool is_integer(std::string inp_str);

bool is_double(std::string inp_str);

bool is_ignore_case_equal(std::string str1, std::string str2);

bool is_single_word(const std::string &inp_str);

size_t n_of_words(const std::string &inp_str);

void skip_n_words(std::stringstream &inp_stream, size_t nskip);

void print_boundary(size_t n);

template <typename T> int sgn(T val);

void time_elapsed( SysClock tick,  SysClock tock);

void debug(int code);
#endif