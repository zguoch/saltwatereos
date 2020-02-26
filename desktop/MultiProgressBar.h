
/**
 * @file MultiProgressBar.h
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Definition of progress bar.
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef MULTIPROGRESSBAR
#define MULTIPROGRESSBAR
#include "stdio.h"
#include <iostream>
#include <vector>
#ifdef _WIN32
        
#else
    // getopt_long only works on MacOS and linux, doesn't work on windows
    #include <sys/ioctl.h>
#endif
#include <stdio.h>
using namespace std;
#include "SWEOSbash.h"
#define COLOR_BAR_PURPLE 0
#define COLOR_BAR_BLUE 1
#define COLOR_BAR_GREEN 2
#define COLOR_BAR_YELLOW 3
#define COLOR_BAR_RED 4

#include <unistd.h>
#define MOVEUP(x) printf("\033[%dA", (x))
// clean screen
#define CLEAR() printf("\033[2J") 
// move up cursor
#define MOVEUP(x) printf("\033[%dA", (x)) 
// move down cursor
#define MOVEDOWN(x) printf("\033[%dB", (x)) 
// move left cursor
#define MOVELEFT(y) printf("\033[%dD", (y)) 
// move right cursor
#define MOVERIGHT(y) printf("\033[%dC",(y)) 
// locate cursor
#define MOVETO(x,y) printf("\033[%d;%dH", (x), (y)) 
// reset cursor
#define RESET_CURSOR() printf("\033[H") 
// hide cursor
#define HIDE_CURSOR() printf("\033[?25l") 
// show cursor
#define SHOW_CURSOR() printf("\033[?25h") 
#define HIGHT_LIGHT() printf("\033[7m")
#define UN_HIGHT_LIGHT() printf("\033[27m")

class MultiProgressBar
{
    public:
        MultiProgressBar(vector<double> xmin,vector<double> xmax,vector<string> title);
        MultiProgressBar(double total,int color=0);
        ~MultiProgressBar();
        void Update(double current_pos=-1);
        void Update(vector<double>current_pos);
    protected:
        vector<string> m_bar_str;
        int m_length_bar;
        char m_bar_char_left;
        char m_bar_char_right;
        vector<double> m_percent;
        vector<string>m_title;
        vector<double> m_total;
        vector<double> m_current_index;
        int m_bar_number;
        vector<double> m_left,m_right;
        int m_maxLength_title;
    private:
        double m_factor;
        vector<string>m_colors;
        void init_colors();
        int m_defaultcolor;
};

#endif //MULTIPROGRESSBAR

