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
    /**
     * @brief Table 6.2 of reference \cite wagner2002iapws.
     * 
     */
    struct Table62
    {
        static const unsigned int numCoeff = 56;
        double c[numCoeff], d[numCoeff], t[numCoeff], n[numCoeff];
        double alpha[numCoeff], beta[numCoeff], gamma[numCoeff], epsilon[numCoeff], a[numCoeff], b[numCoeff], A[numCoeff], B[numCoeff], C[numCoeff], D[numCoeff];
    };
    
    // ============= Constants of H2O ========================================================================
    // double const P_Critic_MPa = 22.064; /**< Pressure at critical point of H2O, [MPa]. IAPWS-95 */
    double const T_Critic_K = 647.096; /**< Temperature at critical point of H2O, [K]. IAPWS-95 */ 
    double const P_Critic = 220.64; /**< Pressure at critical point of H2O, [bar]. IAPWS-95 */
    double const T_Critic = T_Critic_K - Kelvin; /**< Temperature at critical point of H2O, [\f$^{\circ}\text{C} \f$]. IAPWS-95 */ 
    double const Rho_Critic = 322; /**< Critical density, [\f$ kg \ m^{-3} \f$]. IAPWS-95 */ 

    double const MolarMass = 0.018015; /**< molar mass of water  [kg/mol] */

    double const T_Triple_K = 273.16; /**< [K] */
    double const T_Triple = 273.16 - Kelvin; /**< [\f$^{\circ}\text{C} \f$] */
    double const Rho_Triple_liquid = 999.793; /**< Liquid density at triple point, [\f$ kg \ m^{-3} \f$]. IAPWS-95 */ 
    double const Rho_Triple_vapor = 0.00485458; /**< Vapor density at triple point, [\f$ kg \ m^{-3} \f$]. IAPWS-95 */ 

    double const R_const = 0.46151805; /**< Specific gas constant [\f$ kJ kg^{-1} K^{-1} \f$] */
    // =======================================================================================================

    /**
     * @brief Water EOS and properties implementation. 
     * 
     * This implementation is based on public libraries, e.g. <a href="http://freesteam.sourceforge.net">Freesteam 2.0</a>,  <a href="http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas/PROST%20Properties%20of%20Water%20and%20Steam.htm">PROST</a>, and published references \cite huber2009new, \cite wagner2002iapws, \cite Klyukin2020. 
     * 
     * The formulation is valid in the entire stable fluid region from the melting curve to 1273 K at pressures to 1000 MPa. It extrapolates in a physically reasonable way outside this region. See also <a href="http://www.iapws.org/relguide/IAPWS-95.html">IAPWS95</a>.
     * 
     */
    class cH2O : public Fluid::cFluid
    {
    private:
        Table62 m_Table62;
        void LoadTable62(Table62& m_Table62);
        bool m_isHighAccuracy;
    public:
        cH2O(/* args */);
        ~cH2O();
        /**
         * @brief Temperature-pressure relations on boiling curve of water. See equation (2.5) of reference \cite wagner2002iapws.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Pressure [bar]
         */
        virtual double P_Boiling(double T);
        /**
         * @brief Calculate pressure based saturated vapor density #Rho_Vapor_Saturated  and relations between pressure, density and temperautre #Pressure_T_Rho.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Pressure [bar]
         */
        double BoilingCurve(double T);
        /**
         * @brief Saturated liquid density. See equation (2.6) of reference \cite wagner2002iapws. 
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Density [\f$ kg\ m^{-3} \f$]
         */
        double Rho_Liquid_Saturated(double T);
        /**
         * @brief Saturated vapor density. See equation (2.7) of reference \cite wagner2002iapws. 
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Density [\f$ kg\ m^{-3} \f$]
         */
        double Rho_Vapor_Saturated(double T);
        /**
         * @brief Water density as a function of T and P.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Density [\f$ kg\ m^{-3} \f$]
         * @return double Density [\f$ kg\ m^{-3} \f$]
         */
        double Rho(double T, double P);
        /**
         * @brief Derivative of the residual part \f$ \phi^r \f$ of the dimensionless Helmholtz free energy. See Table 6.5 of the reference \cite wagner2002iapws.
         * 
         * 
         * @param delta \f$ \delta = \frac{\rho}{\rho_c} \f$
         * * @param tau \f$ \tau = \frac{T_c}{T} \f$
         * @return double 
         */
        double Phi_r_delta(double delta, double tau);
        /**
         * @brief Calculate pressure given temperature and density. See equation in Table 6.3 of reference \cite wagner2002iapws.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param Rho Density [\f$ kg\ m^{-3} \f$]
         * @return double Pressure [bar]
         */
        double Pressure_T_Rho(double T, double Rho);
    public:
        
    };
    
}

#endif