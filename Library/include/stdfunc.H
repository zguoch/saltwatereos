#ifndef STDFUNC_H
#define STDFUNC_H
#include "stdio.h"
#include<map>  
#ifdef _WIN32
        #include "windows.h"
        #define BLACK			0
        #define BLUE			1
        #define GREEN			2
        #define CYAN			3
        #define RED				4
        #define MAGENTA			5
        #define BROWN			6
        #define LIGHTGRAY		7
        #define DARKGRAY		8
        #define LIGHTBLUE		9
        #define LIGHTGREEN		10
        #define LIGHTCYAN		11
        #define LIGHTRED		12
        #define LIGHTMAGENTA	13
        #define YELLOW			14
        #define WHITE			15
        static HANDLE   m_hConsole=GetStdHandle(STD_OUTPUT_HANDLE);
        static WORD     m_currentConsoleAttr;
        static CONSOLE_SCREEN_BUFFER_INFO csbi;
        #define COLOR_PURPLE ""
        #define COLOR_RED "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<""
        #define COLOR_GREEN "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (GREEN & 0x0F) );cout<<""
        #define COLOR_YELLOW "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<""
        #define COLOR_BLUE "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (BLUE & 0x0F) );cout<<""
        #define COLOR_DEFAULT "";SetConsoleTextAttribute(m_hConsole, m_currentConsoleAttr );cout<<""
        #define ERROR_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (RED & 0x0F) );cout<<"Error: "<<COLOR_DEFAULT
        #define WARN_COUT "";SetConsoleTextAttribute(m_hConsole, ((BLACK & 0x0F) << 4) + (YELLOW & 0x0F) );cout<<"Warning: "<<COLOR_DEFAULT
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
    //-----------------------------------------------------------------------------
#endif