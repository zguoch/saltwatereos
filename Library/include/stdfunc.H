#ifndef STDFUNC_H
#define STDFUNC_H
#include "stdio.h"
#include <iostream>
#include<map>
#include <sys/stat.h>
#include <cstring>
#include <sstream>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

 #ifdef _WIN32
    // define color, this seems only work on MacOS and linux, doesn't work on windows
    #define ERROR_COUT ""
    #define WARN_COUT ""
    #define COLOR_PURPLE ""
    #define COLOR_RED ""
    #define COLOR_GREEN ""
    #define COLOR_YELLOW ""
    #define COLOR_BLUE ""
    #define COLOR_DEFAULT ""
#else
    // define color, this seems only work on MacOS and linux, doesn't work on windows
    #define ERROR_COUT "["<<"\033[31mError: "<<"\033[0m] "
    #define WARN_COUT "["<<"\033[33mWarning: "<<"\033[0m] "
    #define COLOR_PURPLE "\033[35m"
    #define COLOR_RED "\033[31m"
    #define COLOR_GREEN "\033[32m"
    #define COLOR_YELLOW "\033[33m"
    #define COLOR_BLUE "\033[34m"
    #define COLOR_DEFAULT "\033[0m"
#endif

#define STATUS(info) std::cout<<"--  "<<COLOR_GREEN<<info<<COLOR_DEFAULT<<std::endl;
#define STATUS_color(info, color) std::cout<<"--  "<<color<<info<<COLOR_DEFAULT<<std::endl;
#define STATUS_time(info, time_taken) std::cout<<"--  "<<COLOR_GREEN<<info<<", time: "<<(double)(time_taken)/CLOCKS_PER_SEC<<" s"<<COLOR_DEFAULT<<std::endl;
#define STATUS_system_time(info, duration) std::cout<<"--  "<<COLOR_GREEN<<info<<", time: "<< double(duration) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den<<" s"<<COLOR_DEFAULT<<std::endl;
#define ERROR(info) {std::cout<<"--  "<<COLOR_RED<<info<<COLOR_DEFAULT<<std::endl; exit(0);}
#define WARNING(info) {std::cout<<"--  "<<COLOR_YELLOW<<info<<COLOR_DEFAULT<<std::endl;}
#define ASSERT(expression, info) {if(!(expression))std::cout<<"--  "<<COLOR_RED<<info<<COLOR_DEFAULT<<std::endl; assert(expression);}
#define WAIT(where) {std::cout<<"Waiting "<<where<< ". Enter to continue..." << std::endl; std::string dummy; std::getline(std::cin, dummy);}
#endif