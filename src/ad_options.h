#ifndef __ad_options_h__
#define __ad_options_h__

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <algorithm>

char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string & option);
bool read_memory_options(int argc, char** argv, const bool grad_calc);
int read_sub_option2(int argc, char** argv, int cmp_regime);

#endif