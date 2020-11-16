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
//[K]
#define TMIN 273.15  
#define TMAX 1273.15 
/**< [Pa] */
#define PMIN 1e5    
#define PMAX 1000E5 
// ============================================
namespace Fluid
{
    
    /**
     * @brief Base c++ class of fluid Equation of State.
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