#ifndef WATER_H
#define WATER_H
#include<iostream>
#include<fstream>
#include<string> 
#include<cmath>
#include<vector>
using namespace std;
#include "Fluid.H"

namespace H2O
{
    // ============= Constants of H2O ========================================================================
    double const P_Critic = 220.64; /**< Pressure at critical point of H2O, [bar]. IAPWS-95 */
    double const T_Critic = 373.946; /**< Temperature at critical point of H2O, [\f$^{\circ}C \f$]. IAPWS-95 */ 
    double const MolarMass = 0.018015; /**< molar mass of water  [kg/mol] */
    // =======================================================================================================

    /**
     * @brief Water EOS and properties implementation. 
     * 
     * This implementation is based on public libraries, e.g. <a href="http://freesteam.sourceforge.net">Freesteam 2.0</a>,  <a href="http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas/PROST%20Properties%20of%20Water%20and%20Steam.htm">PROST</a>, and published references \cite huber2009new, \cite wagner2002iapws, \cite Klyukin2020. 
     * 
     */
    class cH2O : public Fluid::cFluid
    {
    private:
        /* data */
    public:
        cH2O(/* args */);
        ~cH2O();
        /**
         * @brief Temperature-pressure relations on boiling curve of water. See equation (2.5) of reference \cite wagner2002iapws.
         * 
         * @param T Temperature [\f$ ^{\circ}C \f$]
         * @return double Pressure [bar]
         */
        virtual double P_Boiling(double T);
    public:
        
    };
    
}

#endif