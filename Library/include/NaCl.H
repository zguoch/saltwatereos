#ifndef NaCl_H
#define NaCl_H
#include<iostream>
#include<fstream>
#include<string> 
#include<cmath>
#include<vector>
using namespace std;
#include "Fluid.H"

namespace NaCl
{
    // ============= Constants of NaCl =======================================
    double const MolarMass = 0.058443; /**< molar mass of water  [kg/mol] */
    double const T_Triple = 800.7; /**< Temperature at triple point of NaCl, [\f$^{\circ}C \f$] */ 
    double const P_Triple = 0.0005; /**< Pressure at triple point of NaCl, [Pa] */
    // =======================================================================

    /**
     * @brief Water EOS and properties function.
     * 
     */
    class cNaCl : public Fluid::cFluid
    {
    private:
        /* data */
    public:
        cNaCl(/* args */);
        ~cNaCl();
        /**
         * @brief T-P relations of boling curve of NaCl.
         * Equation (3) of reference \cite Driesner2007Part1. 
         * 
         * \f$ log_{10}(P_{NaCl, liquid}) = log_{10}(P_{triple, liquid}) + b_{boil}\left( \frac{1}{T_{triple, NaCl} + 273.15} - \frac{1}{T + 273.15} \right) , b_{boil} = 0.941812\times 10^4 \f$
         * 
         * @param T Temperature, [\f$ ^{\circ}C\f$]
         * @return double Sublimation pressure of NaCl, [bar]
         * 
         * \image html Driesner_Heinrich_Fig4.png "Vapor pressure curves (sublimation and boiling curves) of NaCl." width=25%.
         * Vapor pressure curves (sublimation and boiling curves) of NaCl, figure 4 of reference \cite Driesner2007Part1. 
         * 
         * \image html HaliteSublimationBoilingCurves.svg "Vapor pressure curves (sublimation and boiling curves) of NaCl." width=25%.
         * see also #P_Sublimation 
         */
        virtual double P_Boiling(double T);
        
        // virtual double P_Boiling(double T);
        /**
         * @brief Melting curve of NaCl. 
         * Equation (1) of reference \cite Driesner2007Part1. 
         * 
         * \f$ T_{hm} = T_{triple, NaCl} + a(P - P_{triple, NaCl}), a=2.47260\times 10^{-2} \f$
         * @param p Pressure, [Pa]
         * @return double Melting temperature, [\f$ ^{\circ}C\f$]
         * 
         * \image html Driesner_Heinrich_Fig3.png "Halite melting curve." width=25%.
         * Halite melting curve, figure 3 of reference \cite Driesner2007Part1. 
         * 
         * \image html HaliteMeltingCurve.svg "Halite melting curve calculated by swEOS." width=25%. 
         */
        inline double T_Melting(double P){ return T_Triple + 2.4726e-2*(P - P_Triple);};
        /**
         * @brief Halite vapor pressure as a function of temperature T, called halite sublimation curve as well.
         * Equation (2) of reference \cite Driesner2007Part1. 
         * 
         * \f$ log_{10}(P_{NaCl, halite}) = log_{10}(P_{triple, NaCl}) + b_{subl}\left( \frac{1}{T_{triple, NaCl} + 273.15} - \frac{1}{T + 273.15} \right) , b_{subl} = 1.18061\times 10^4 \f$
         * 
         * @param P Temperature, [\f$ ^{\circ}C\f$]
         * @return double Sublimation pressure of NaCl, [bar]
         * 
         * \image html Driesner_Heinrich_Fig4.png "Vapor pressure curves (sublimation and boiling curves) of NaCl." width=25%.
         * Vapor pressure curves (sublimation and boiling curves) of NaCl, figure 4 of reference \cite Driesner2007Part1. 
         * 
         * \image html HaliteSublimationBoilingCurves.svg "Vapor pressure curves (sublimation and boiling curves) of NaCl." width=25%. 
         * see also #P_Boiling
         */
        inline double P_Sublimation(double T){return pow(10, (log10(P_Triple) + 11806.1 * (1.0 / (T_Triple + 273.15) - 1.0 / (T + 273.15))));};
     
    public:
        
    };
    
}

#endif