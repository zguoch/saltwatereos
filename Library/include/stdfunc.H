#ifndef STDFUNC_H
#define STDFUNC_H
#include "stdio.h"
#include<map>  

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

#endif