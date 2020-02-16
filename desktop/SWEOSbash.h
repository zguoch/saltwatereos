#ifndef SWEOSBASH_H
#define SWEOSBASH_H

#include <unistd.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <sys/ioctl.h>
#include <vector>
using namespace std;

#include "H2ONaCl.H"



namespace SWEOSbash
{
    #ifdef _WIN32
        #include "windows.h"
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

    void bash_run(int argc, char** argv);
    // calculation mode: 0d, 1d, 2d, 3d
    #define CALCULATION_MODE_SINGLEPOINT 0
    #define CALCULATION_MODE_ONEDIMENSION 1
    #define CALCULATION_MODE_TWODIMENSION 2
    #define CALCULATION_MODE_THREEDIMENSION 3
    
    // variable selection
    #define VARIABLE_SELECTION_PTX 0
    #define VARIABLE_SELECTION_PHX 1
    #define VARIABLE_SELECTION_T 2
    #define VARIABLE_SELECTION_P 3
    #define VARIABLE_SELECTION_X 4
    #define VARIABLE_SELECTION_H 5
    #define VARIABLE_SELECTION_PT 6
    #define VARIABLE_SELECTION_PX 7
    #define VARIABLE_SELECTION_TX 8
    #define VARIABLE_SELECTION_PH 9
    #define VARIABLE_SELECTION_HX 10

    class cSWEOSarg
    {
    private:
        bool m_haveD, m_haveV, m_haveP, m_haveT, m_haveX, m_haveH, m_haveR;
        int m_valueD;
        string m_valueV;
        double m_valueT, m_valueP, m_valueX, m_valueH;
        // min/delta/max, order coresponding to -V parameter, 
        //e.g. -VPT, m_valueR1 for pressure, m_valueR2 for temperature
        // double m_valueR1[3], m_valueR2[3], m_valueR3[3];
        double m_valueR[3][3];
        vector<string> m_varueR_str;
    public:
        cSWEOSarg(/* args */);
        ~cSWEOSarg();
    public:
        int m_CalculationMode;
        int m_VariableSelection;
    public:
        bool ParseAndCheck(int argc, char** argv);
    private:
        vector<string> string_split(string s, string delimiter);
    };

    SWEOS::PROP_H2ONaCl calculateSinglePoint_PTX(double P, double T, double X, bool isCout=true);
    SWEOS::PROP_H2ONaCl calculateSinglePoint_PHX(double P, double H, double X, bool isCout=true);
    static void StartText()
    {
        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<COLOR_YELLOW;       //print text in yellow color
        cout << "***************************************************\n";
        cout << "*                 program swEOS                   *\n";
        cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
        cout << "*                                                 *\n";
        cout << "*  Equation of state of salt water                *\n";
        cout << "*  - Independent variables: PTX, PHX              *\n";
        cout << "*  - Properties: density, enthalpy, viscosity     *\n";
        cout << "*  - saturation, salinity, phase diagram          *\n";
        cout << "*  unit:                                          *\n";
        cout << "*      temperature-°C,        pressure-bar        *\n";
        cout << "*      salinity-wt. % NaCl,   density-kg/m3       *\n";
        cout << "*      enthalpy-kJ/kg,        viscosity-Pa s      *\n";
        cout << "*                                                 *\n";
        cout << "* (c) Zhikui Guo, GEOMAR, Feb 15 2020, Kiel       *\n";
        cout << "*                                                 *\n";
        cout << "***************************************************\n";
        cout << "\n";
        cout<<COLOR_DEFAULT;
                                                                                                                                                                                                                        
    }
    static void helpINFO()
    {
        string version="1.0";
        string author="Zhikui Guo";
        string locus="GEOMAR, Germany";
        unsigned int wordWidth=20;
        // time_t now=time(0);
        // char* now_str=ctime(&now);
        string now_str="Feb 15, 2020";

        //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
        cout<<"========================== swEOS ==========================="<<endl;
        cout<<"Analytical continuation of potential field data"<<endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author "<<COLOR_GREEN<<author<<COLOR_DEFAULT<<endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus "<<COLOR_GREEN<<locus<<COLOR_DEFAULT<<endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date "<<COLOR_GREEN<<now_str<<COLOR_DEFAULT<<endl;
        cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version "<<COLOR_GREEN<<version<<COLOR_DEFAULT<<endl;
        cout<<"============================================================"<<endl;
    }
    static void StartText_artASCII()
    {
        // cout<<COLOR_GREEN<<"███████╗ █████╗ ██╗  ████████╗██╗    ██╗ █████╗ ████████╗███████╗██████╗     ███████╗ ██████╗ ███████╗\n"
        // <<"██╔════╝██╔══██╗██║  ╚══██╔══╝██║    ██║██╔══██╗╚══██╔══╝██╔════╝██╔══██╗    ██╔════╝██╔═══██╗██╔════╝\n"
        // <<"███████╗███████║██║     ██║   ██║ █╗ ██║███████║   ██║   █████╗  ██████╔╝    █████╗  ██║   ██║███████╗\n"
        // <<"╚════██║██╔══██║██║     ██║   ██║███╗██║██╔══██║   ██║   ██╔══╝  ██╔══██╗    ██╔══╝  ██║   ██║╚════██║\n"
        // <<"███████║██║  ██║███████╗██║   ╚███╔███╔╝██║  ██║   ██║   ███████╗██║  ██║    ███████╗╚██████╔╝███████║\n"
        // <<"╚══════╝╚═╝  ╚═╝╚══════╝╚═╝    ╚══╝╚══╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚═╝  ╚═╝    ╚══════╝ ╚═════╝ ╚══════╝\n"
        // <<COLOR_DEFAULT<<endl;                                                                                                 

        cout<<COLOR_GREEN<<"                              $$$$$$$$\\  $$$$$$\\   $$$$$$\\  \n"
        <<"                              $$  _____|$$  __$$\\ $$  __$$\\ \n"
        <<" $$$$$$$\\ $$\\  $$\\  $$\\       $$ |      $$ /  $$ |$$ /  \\__|\n"
        <<"$$  _____|$$ | $$ | $$ |      $$$$$\\    $$ |  $$ |\\$$$$$$\\  \n"
        <<"\\$$$$$$\\  $$ | $$ | $$ |      $$  __|   $$ |  $$ | \\____$$\\ \n"
        <<" \\____$$\\ $$ | $$ | $$ |      $$ |      $$ |  $$ |$$\\   $$ |\n"
        <<"$$$$$$$  |\\$$$$$\\$$$$  |      $$$$$$$$\\  $$$$$$  |\\$$$$$$  |\n"
        <<"\\_______/  \\_____\\____/       \\________| \\______/  \\______/ \n"
        <<COLOR_DEFAULT<<endl;                                                                                                          
    }
}
#endif