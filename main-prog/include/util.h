#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <map>

void print_banner();

void debug(std::string msg);

bool has_input(const int inp_num);

bool is_integer(std::string inp_str);

bool is_double(std::string inp_str);

bool is_ignore_case_equal(std::string str1, std::string str2);

bool is_single_word(const std::string &inp_str);

unsigned int n_of_words(const std::string &inp_str);

void skip_n_words(std::stringstream &inp_stream, unsigned int nskip);

void print_boundary(unsigned int n);

template <typename T> int sgn(T val);

#endif