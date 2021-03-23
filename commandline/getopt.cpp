#include "getopt.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <string.h>
#endif

#include <iostream>

char* optarg = NULL;
int optind = 1;

int getopt(int argc, char *const argv[], const char *optstring)
{
    if (optind >= argc) //check if the arguments are parsed at end
    {
        return -1;
    }
    if ((optind >= argc) || (argv[optind][0] != '-') || (argv[optind][0] == 0))
    {
        return -1;
    }
    int opt = argv[optind][1];
    const char *p = strchr(optstring, opt);
    if (p == NULL)
    {
        optind++;
        return '?';
    }
    if (p[1] == ':') //argument with option
    {
        optind++;
        if (optind >= argc)
        {
            return '?';
        }
        optarg = argv[optind];
        optind++;
    }else //argument without option
    {
        optind++;
    }
    
    return opt;
}
