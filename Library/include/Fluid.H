#ifndef FLUID_H
#define FLUID_H
#include<iostream>
#include<fstream>
#include<string> 
#include<cmath>
#include<vector>
using namespace std;

// ============= Constants ====================
#define Kelvin 273.15
// ============================================
namespace Fluid
{
    
    /**
     * @brief Base C++ class of fluid Equation of State. 
     * 
     * EOS and properties functions are designed as pure virtual functions and then implemented in the specific specie classes, e.g. H2O::cH2O, NaCl::cNaCl, H2ONaCl::cH2ONaCl. While some commonly used functions are defined and implemented in the base class directly. 
     * 
     */
    class cFluid
    {
    private:
        /* data */
    public:
        cFluid(/* args */);
        ~cFluid();
        /**
         * @brief Boiling curve.
         * 
         * @param T T Temperature, [\f$ ^{\circ}C\f$]
         * @return double Pressure [bar]
         */
        virtual double P_Boiling(double T) = 0;
    };
    
}

#endif