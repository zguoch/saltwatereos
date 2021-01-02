#ifndef WATER_H
#define WATER_H
#include <iostream>
#include <fstream>
#include <string> 
#include <cmath>
#include <vector>
using namespace std;
// #include "Fluid.H"
// ============= Constants ====================
#define Kelvin 273.15
// ============================================
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
    /**
     * @brief Table 6.1 of \cite wagner2002iapws.
     * 
     */
    struct Table61
    {
        double n0[8]={-8.32044648201, 6.6832105268, 3.00632, 0.012436, 0.97315, 1.2795, 0.96956, 0.24873};
        double gamma0[8]={0, 0, 0, 1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105};
    };
    /**
     * @brief Table 2, 3 and 5 of \cite huber2009new.
     * 
     */
    struct Table235
    {
        double H[4]={1.67752, 2.20462, 0.6366564, -0.241605}; //table 2
        double Hij[6][7]=
        {
            {0.520094,  0.222531,   -0.281378,  0.161913,   -0.0325372,     0,              0},
            {0.0850895, 0.999115,   -0.906851,  0.257399,   0,              0,              0},
            {-1.08374,  1.88797,    -0.772479,  0,          0,              0,              0},
            {-0.289555, 1.26613,    -0.489837,  0,          0.0698452,      0,              -0.00435673},
            {0,         0,          -0.25704,   0,          0,              0.00872102,     0},
            {0,         0.120573,   0,          0,          0,              0,              -0.000593264}
        }; //table 3
        double chi_mu=0.068, q_C=1.0/1.9, q_D=1.0/1.1, nu=0.63, gamma=1.239, xi0=0.13, Gamma0=0.06, T_bar_R=1.5;
    };
    
    // ============= Constants of H2O ========================================================================
    double const PMIN = 1E5; /**< Minimum valid pressure of H2O, [Pa]. IAPWS-95 */ 
    double const PMAX = 10000E5; /**< Minimum valid pressure of H2O, [Pa]. IAPWS-95 */ 
    double const TMIN = 0.1; /**< Minimum valid pressure of H2O, [\f$ ^{\circ}\text{C} \f$]. IAPWS-95 */ 
    double const TMAX = 1000; /**< Minimum valid pressure of H2O, [\f$ ^{\circ}\text{C} \f$]. IAPWS-95 */ 
    // double const P_Critic_MPa = 22.064; /**< Pressure at critical point of H2O, [MPa]. IAPWS-95 */
    double const T_Critic_K = 647.096; /**< Temperature at critical point of H2O, [K]. IAPWS-95 */ 
    double const P_Critic = 220.64; /**< Pressure at critical point of H2O, [bar]. IAPWS-95 */
    double const T_Critic = T_Critic_K - Kelvin; /**< Temperature at critical point of H2O, [\f$^{\circ}\text{C} \f$]. IAPWS-95 */ 
    double const Rho_Critic = 322; /**< Critical density, [\f$ kg \ m^{-3} \f$]. IAPWS-95 */ 

    double const MolarMass = 0.018015; /**< molar mass of water  [kg/mol] */

    double const T_Triple_K = 273.16; /**< [K] */
    double const T_Triple = 273.16 - Kelvin; /**< [\f$^{\circ}\text{C} \f$] */
    double const P_Triple = 611.655E-5; /**< Pressure at triple point of H2O, [bar]. IAPWS-95 */
    double const Rho_Triple_liquid = 999.793; /**< Liquid density at triple point, [\f$ kg \ m^{-3} \f$]. IAPWS-95 */ 
    double const Rho_Triple_vapor = 0.00485458; /**< Vapor density at triple point, [\f$ kg \ m^{-3} \f$]. IAPWS-95 */ 

    enum PhaseRegion {iceI, iceIII, iceV, iceVI, iceVII};
    double const T_K_ice_min[5]={251.165, 251.165, 256.164, 273.31, 355};
    double const T_K_ice_max[5]={T_Triple_K, 256.164, 273.31, 355, 715};

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
    class cH2O// : public Fluid::cFluid
    {
    private:
        Table62 m_Table62;  /*< table 6.2 of wagner(2002) */
        Table61 m_Table61; //table 6.1 of wagner(2002)
        Table235 m_Table235; //table 2 and table 3 of Huber(2009)
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
        double P_Boiling(double T);
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
         * @param P Pressure [\f$ bar \f$]
         * @return double Density [\f$ kg\ m^{-3} \f$]
         */
        double Rho(double T, double P);
        /**
         * @brief See equation 6.6 of \cite wagner2002iapws. 
         * 
         * \f{matrix}
         * \phi^r & = & \sum\limits_{i=1}^{7} n_i \delta^{d_i}\tau^{t_i} + \sum\limits_{i=8}^{51} n_i \delta^{d_i} e^{-\delta^{c_i}}  \tau^{t_i} + \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 1} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} + \sum\limits_{i=55}^{56} n_i \delta \Delta^{b_i}\psi \\
         * \Delta & = & \theta^2 + B_i [\color{red}{(\delta -1)^2}]^{a_i} \\
         * \theta & = & (1-\tau) + A_i [\color{red}{(\delta -1)^2}]^{1/2\beta_i} \\
         * \psi & = & e^{-C_i (\delta -1)^2 - D_i (\tau - 1)^2}
         * \f}
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * \note \f$\color{red}{ (\delta -1)^2} \f$ term has to be kept to make sure it is a positive 
         * 
         * @param delta \f$ \delta =  \rho / \rho_c\f$
         * @param tau \f$ \tau = T_c / T\f$, unit of \f$ T, T_c \f$ is K.
         * @return double 
         */
        double Phi_r(double delta, double tau);
        /**
         * @brief Derivative of the residual part \f$ \phi^r \f$ of the dimensionless Helmholtz free energy. See Table 6.5 of the reference \cite wagner2002iapws.
         * 
         * \f$ \phi^{r}_{\delta} = \left[ \frac{\partial \phi^r}{\partial \delta} \right]_{\tau} \f$.  
         * \f$ \phi^r\f$  is defined in equation 6.6 of \cite wagner2002iapws, where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * \f{equation}
         * \phi^r_{\delta} = part_1 + part_2 + part_3 + part_4
         * \f} 
         * 
         * \f{equation}
         * \begin{matrix}
         * part_1 & = & \sum\limits_{i=1}^{7} n_i d_i \delta^{d_i -1} \tau^{t_i} \\
         * part_2 & = & \sum\limits_{i=8}^{51} n_i \delta^{d_i - 1}  e^{-\delta^{c_i}}  \tau^{t_i} (d_i -c_i \delta^{c_i}) \\
         * part_3 & = & \sum\limits_{i=52}^{54} n_i \delta^{d_i - 1} \tau^{t_i} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[ d_i -2\delta \alpha_i  (\delta - \epsilon_i) \right] \\
         * part_4 & = & \sum\limits_{i=55}^{56} n_i \left[ \Delta^{b_i}\psi +  \delta \left( \psi\color{blue}{\frac{\partial \Delta^{bi}}{\partial \delta}} + \Delta^{b_i}\color{blue}{\frac{\partial \psi}{\partial \delta}} \right) \right] \\
         * \color{blue}{\frac{\partial \Delta^{bi}}{\partial \delta}} & = & 2b_i\Delta^{b_i -1} (\delta -1) \left[ \theta \frac{A_i}{\beta_i} [ \color{red}{(\delta -1)^2} ]^{1/2\beta_i -1} + B_ia_i[\color{red}{ (\delta -1)^2 } ]^{a_i-1} \right] \\
         * \color{blue}{\frac{\partial \psi}{\partial \delta}} & = & -2C_i(\delta -1) \psi 
         * \end{matrix}
         * \f} 
         * \note \f$\color{red}{ (\delta -1)^2} \f$ term has to be kept to make sure it is a positive value.
         * @param delta \f$ \delta = \frac{\rho}{\rho_c} \f$
         * @param tau \f$ \tau = \frac{T_c}{T} \f$
         * @return double 
         */
        double Phi_r_delta(double delta, double tau);
        /**
         * @brief Derivative of the residual part \f$ \phi^r \f$ of the dimensionless Helmholtz free energy. See Table 6.5 of the reference \cite wagner2002iapws.
         * 
         * \f$ \phi^{r}_{\delta\delta} = \left[ \frac{\partial^2 \phi^r}{\partial \delta^2} \right]_{\tau} \f$.  
         * \f$ \phi^r\f$  is defined in equation 6.6 of \cite wagner2002iapws, where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * \f{equation}
         * \phi^r_{\delta\delta} = part_1 + part_2 + part_3 + part_4
         * \f} 
         * 
         * \f{equation}
         * \begin{matrix}
         * part_1 & = & \sum\limits_{i=1}^{7} n_i d_i (d_i - 1) \delta^{d_i -2} \tau^{t_i} \\
         * part_2 & = & \sum\limits_{i=8}^{51} n_i \delta^{d_i - 2}  e^{-\delta^{c_i}}  \tau^{t_i}\left[ (d_i -c_i \delta^{c_i})(d_i -1 -c_i \delta^{c_i}) - c_i^2\delta^{c_i} \right] \\
         * part_3 & = & \sum\limits_{i=52}^{54} n_i \delta^{d_i - 2} \tau^{t_i} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[  \color{red}{[d_i - 1 -2\delta\alpha_i(\delta-\epsilon_i)]}\color{green}{[d_i -2\delta\alpha_i(\delta-\epsilon_i)]} -\color{blue}{\delta (4\delta\alpha_i + 2\alpha_i \epsilon_i)} \right] \\
         * part_4 & = & \sum\limits_{i=55}^{56} n_i \left[ \color{blue}{\frac{\partial \Delta^{b_i}}{\partial \delta}}\psi + \Delta^{b_i}\frac{\partial \psi}{\partial \delta} +  \left( \psi\color{blue}{\frac{\partial \Delta^{bi}}{\partial \delta}} + \Delta^{b_i}\color{blue}{\frac{\partial \psi}{\partial \delta}} \right) +\delta \left(\frac{\partial \psi}{\partial \delta}\color{blue}{\frac{\partial \Delta^{bi}}{\partial \delta}}  +  \psi\color{blue}{\frac{\partial^2 \Delta^{bi}}{\partial \delta^2}} + \frac{\partial \Delta^{b_i}}{\partial \delta}\color{blue}{\frac{\partial \psi}{\partial \delta}} + \Delta^{b_i}\color{blue}{\frac{\partial^2 \psi}{\partial \delta^2}} \right) \right] \\
         * \color{blue}{\frac{\partial \Delta^{bi}}{\partial \delta}} & = & b_i\Delta^{b_i -1} \frac{\partial \Delta}{\partial \delta} \\
         * \frac{\partial^2 \Delta^{b_i}}{\partial \delta^2} & = & b_i\Delta^{b_i - 2} \left[ (b_i -1)(\frac{\partial \Delta}{\partial \delta})^2 + \Delta \frac{\partial^2 \Delta}{\partial \delta^2} \right] \\
         * \color{blue}{\frac{\partial \Delta}{\partial \delta}} & = & 2(\delta -1) \left[ \theta \frac{A_i}{\beta_i} [ \color{red}{(\delta -1)^2} ]^{1/2\beta_i -1} + B_ia_i[\color{red}{ (\delta -1)^2 } ]^{a_i-1} \right] \\
         * \color{blue}{\frac{\partial \psi}{\partial \delta}} & = & -2C_i(\delta -1) \psi 
         * \end{matrix}
         * \f} 
         * \note \f$\color{red}{ (\delta -1)^2} \f$ term has to be kept to make sure it is a positive value.
         * @param delta \f$ \delta = \frac{\rho}{\rho_c} \f$
         * @param tau \f$ \tau = \frac{T_c}{T} \f$
         * @return double 
         */
        double Phi_r_deltadelta(double delta, double tau);
        /**
         * @brief See equation 6.5 of \cite wagner2002iapws. 
         * 
         * \f{equation}
         * \phi^o = ln\delta + n_1^o + n_2^o\tau + n_3^o ln\tau + \sum\limits_{i=4}^8 n_i^o ln[1-e^{-\gamma_i^o \tau}]
         * \f}
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * @param delta \f$ \delta =  \rho / \rho_c\f$
         * @param tau \f$ \tau = T_c / T\f$, unit of \f$ T, T_c \f$ is K.
         * @return double 
         */
        double Phi_o(double delta, double tau);
        /**
         * @brief See Table 6.3 of \cite wagner2002iapws. 
         * 
         * \f$ \phi^{o}_{\tau} = \left[ \frac{\partial \phi^0}{\partial \tau} \right]_{\delta} \f$.  
         * \f$ \phi^o\f$  is defined in equation 6.5 of \cite wagner2002iapws.
         * \f{equation}
         * \phi^o_{\tau} = n^o_2 + \frac{n^o_3}{\tau} + \sum\limits_{i=4}^8 n^o_i \gamma^o_i \left( \frac{1}{1-e^{-\gamma^o_i \tau}} - 1 \right)
         * \f} 
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * @param delta \f$ \delta =  \rho / \rho_c\f$
         * @param tau \f$ \tau = T_c / T\f$, unit of \f$ T, T_c \f$ is K.
         * @return double 
         */
        double Phi_o_tau(double delta, double tau);
        /**
         * @brief See Table 6.3 of \cite wagner2002iapws. 
         * 
         * \f$ \phi^{r}_{\tau\tau} = \left[ \frac{\partial \phi^r}{\partial \tau} \right]_{\delta} \f$.  
         * \f$ \phi^r\f$  is defined in equation 6.6 of \cite wagner2002iapws.
         * \f{equation}
         * \phi^r_{\tau} = part_1 + part_2 + part_3 + part_4
         * \f} 
         * 
         * \f{equation}
         * \begin{matrix}
         * part_1 = & \sum\limits_{i=1}^{7} n_i \delta^{d_i}t_i \tau^{t_i -1} \\
         * part_2 = & \sum\limits_{i=8}^{51} n_i \delta^{d_i} t_i e^{-\delta^{c_i}}  \tau^{t_i-1} \\
         * part_3 = & \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 1} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[ t_i -2\tau \beta_i  (\tau - \gamma_i) \right] \\
         * part_4 = & \sum\limits_{i=55}^{56} n_i \delta \left( \psi\color{blue}{\frac{\partial \Delta^{bi}}{\partial \tau}} + \Delta^{b_i}\color{blue}{\frac{\partial \psi}{\partial \tau}} \right) \\
         * \color{blue}{\frac{\partial \Delta^{bi}}{\partial \tau}} = & -2b_i\theta\Delta^{b_i - 1} \\
         * \color{blue}{\frac{\partial \psi}{\partial \tau}} = & -2D_i(\tau -1) \psi
         * \end{matrix}
         * \f} 
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * @param delta \f$ \delta =  \rho / \rho_c\f$
         * @param tau \f$ \tau = T_c / T\f$, unit of \f$ T, T_c \f$ is K.
         * @return double 
         */
        double Phi_r_tau(double delta, double tau);
        /**
         * @brief See Table 6.3 of \cite wagner2002iapws. 
         * 
         * \f$ \phi^{o}_{\tau\tau} = \left[ \frac{\partial^2 \phi^o}{\partial \tau ^2} \right]_{\delta} \f$.  
         * \f$ \phi^o\f$  is defined in equation 6.5 of \cite wagner2002iapws.
         * \f{equation}
         * \phi^o_{\tau\tau} = - \frac{n^o_3}{\tau^2} - \sum\limits_{i=4}^8  \frac{n^o_i \gamma^{o2}_i e^{-\gamma^o_i \tau}}{1-e^{-\gamma^o_i \tau}}  
         * \f} 
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * @param delta \f$ \delta =  \rho / \rho_c\f$
         * @param tau \f$ \tau = T_c / T\f$, unit of \f$ T, T_c \f$ is K.
         * @return double 
         */
        double Phi_o_tautau(double delta, double tau);
        /**
         * @brief See Table 6.3 of \cite wagner2002iapws. 
         * 
         * \f$ \phi^{r}_{\tau\tau} = \left[ \frac{\partial^2 \phi^r}{\partial \tau^2} \right]_{\delta} \f$.  
         * \f$ \phi^r\f$  is defined in equation 6.6 of \cite wagner2002iapws.
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * @brief See Table 6.3 of \cite wagner2002iapws. 
         * 
         * \f$ \phi^{r}_{\tau\tau} = \left[ \frac{\partial \phi^r}{\partial \tau} \right]_{\delta} \f$.  
         * \f$ \phi^r\f$  is defined in equation 6.6 of \cite wagner2002iapws.
         * \f{equation}
         * \phi^r_{\tau\tau} = part_1 + part_2 + part_3 + part_4
         * \f} 
         * 
         * \f{equation}
         * \begin{matrix}
         * part_1 = & \sum\limits_{i=1}^{7} n_i \delta^{d_i}t_i(t_i -1) \tau^{t_i -2} \\
         * part_2 = & \sum\limits_{i=8}^{51} n_i \delta^{d_i} t_i(t_i - 1) e^{-\delta^{c_i}}  \tau^{t_i-2} \\
         * part_3 = & \left ( \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 1} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[ t_i -2\tau \beta_i  (\tau - \gamma_i) \right] \right)_{\tau} \\
         * = & \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 2} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[ t_i -1 -2\tau \beta_i  (\tau - \gamma_i) \right] \left[ t_i -2\tau \beta_i  (\tau - \gamma_i) \right]  + n_i \delta^{d_i} \tau^{t_i - 1} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} (2\beta_i\gamma_i - 4\tau\beta_i) \\
         * = & \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 2} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[ \left[ t_i -2\tau \beta_i  (\tau - \gamma_i) \right]^2 -(t_i -2\tau \beta_i  (\tau - \gamma_i)) + \tau(2\beta_i\gamma_i - 4\tau\beta_i) \right] \\
         * = & \sum\limits_{i=52}^{54} n_i \delta^{d_i} \tau^{t_i - 2} e^{-\alpha_i (\delta -\epsilon_i)^2 - \beta_i(\tau - \gamma_i)^2} \left[ \left[ t_i -2\tau \beta_i  (\tau - \gamma_i) \right]^2 -t_i - 2\tau^2\beta_i\right] \\
         * part_4 = & \left ( \sum\limits_{i=55}^{56} n_i \delta \left( \psi\color{blue}{\frac{\partial \Delta^{b_i}}{\partial \tau}} + \Delta^{b_i}\color{blue}{\frac{\partial \psi}{\partial \tau}} \right) \right)_{\tau} \\
         * = & \sum\limits_{i=55}^{56} n_i \delta \left( \color{blue}{\frac{\partial \psi}{\partial \tau}}\color{blue}{\frac{\partial \Delta^{b_i}}{\partial \tau}} + \psi\color{blue}{\frac{\partial^2 \Delta^{b_i}}{\partial \tau^2}}  + \color{blue}{\frac{\partial \Delta^{b_i}}{\partial \tau}}\color{blue}{\frac{\partial \psi}{\partial \tau}} + \Delta^{b_i}\color{blue}{\frac{\partial^2 \psi}{\partial \tau^2}} \right) \\
         * \color{blue}{\frac{\partial \Delta^{bi}}{\partial \tau}} = & -2b_i\theta\Delta^{b_i - 1} \\
         * \color{blue}{\frac{\partial \psi}{\partial \tau}} = & -2D_i(\tau -1) \psi \\
         * \color{blue}{\frac{\partial^2 \Delta^{bi}}{\partial \tau^2}} = & 2b_i\Delta^{b_i-2}\left[ \Delta + 2\theta^2 (b_i -1) \right] \\
         * \color{blue}{\frac{\partial^2 \psi}{\partial \tau^2}} = & 2D_i \psi \left[2D_i(\tau -1)^2 - 1 \right]
         * \end{matrix}
         * \f} 
         * Where \f$ \delta, \tau \f$ are defined in equation 6.4 of \cite wagner2002iapws.
         * 
         * @param delta \f$ \delta =  \rho / \rho_c\f$
         * @param tau \f$ \tau = T_c / T\f$, unit of \f$ T, T_c \f$ is K.
         * @return double 
         */
        double Phi_r_tautau(double delta, double tau);
        double Phi_r_deltatau(double delta, double tau);
        /**
         * @brief Calculate pressure given temperature and density. See equation in Table 6.3 of reference \cite wagner2002iapws.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param Rho Density [\f$ kg\ m^{-3} \f$]
         * @return double Pressure [bar]
         */
        double Pressure_T_Rho(double T, double Rho);
        /**
         * @brief Sublimation-Pressure. \b Solid-Vapor phase boundary. The valid temperature region is \f$[T_{min}, T_{triple}] \f$. See equation (2.21) of reference \cite wagner2002iapws.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Pressure [bar]
         */
        double SublimationCurve(double T);
        /**
         * @brief \b Liquid-Vapor phase boundary. Calculate pressure based saturated vapor density #Rho_Vapor_Saturated  and relations between pressure, density and temperautre #Pressure_T_Rho. The valid temperature region is \f$[T_{triple}, T_{critic}] \f$.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Pressure [bar]
         */
        double BoilingCurve(double T);
        /**
         * @brief Melting-Pressure. \b Solid-Liquid phase boundary. See equation (2.16) of reference \cite wagner2002iapws.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Pressure [bar]
         */
        double MeltingCurve(double T, bool isIceI=false);
        /**
         * @brief Specific enthalpy as a function of temperature and density.
         * See Table 6.3 of \cite wagner2002iapws.
         * 
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param Rho Density [\f$ kg/m^3 \f$]
         * @return Specific enthalpy [\f$ kJ/kg \f$]
         */
        double SpecificEnthalpy_T_Rho(double T, double Rho);
        /**
         * @brief 
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [\f$ bar \f$]
         * @return Specific enthalpy [\f$ kJ/kg \f$] 
         */
        double SpecificEnthalpy(double T, double P);
        /**
         * @brief Isochoric heat capacity as a function of temperature and density.
         * See Table 6.3 of \cite wagner2002iapws.
         * 
         * \f{equation}
         * c_v(T,p)=c_v(T,\rho)=c_v(\delta, \tau) = -R\tau^2 (\phi_{\tau\tau}^o + \phi_{\tau\tau}^r)
         * \f}
         * 
         * \f$\phi_{\tau\tau}^o\f$ is calculated from H2O::cH2O::Phi_o_tautau, and \f$\phi_{\tau\tau}^r\f$ is calculated from H2O::cH2O::Phi_r_tautau
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param Rho Density [\f$ kg/m^3 \f$]
         * @return Isochoric heat capacity [\f$ kJ/kg/K \f$]
         */
        double Cv_T_Rho(double T, double Rho);
        /**
         * @brief Isochoric heat capacity
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [\f$ bar \f$]
         * @return Isochoric heat capacity [\f$ kJ/kg/K \f$] 
         */
        double Cv(double T, double P);
        /**
         * @brief Isobaric heat capacity
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param Rho Density [\f$ kg/m^3 \f$]
         * @return Isobaric heat capacity [\f$ kJ/kg/K \f$] 
         */
        double Cp_T_Rho(double T, double Rho);
        /**
         * @brief Isobaric heat capacity
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [\f$ bar \f$]
         * @return Isobaric heat capacity [\f$ kJ/kg/K \f$] 
         */
        double Cp(double T, double P);
        /**
         * @brief Dynamic viscosity.
         * 
         * See equation 2, 11, 16, 20, 26, 33 and table 2, 3, 5  of \cite huber2009new.
         * 
         * \f{equation}
         * \bar\mu = \bar{\mu_0}(\bar{T}) \times \bar{\mu_1}(\bar{T}, \bar{\rho}) \times \bar{\mu_2}(\bar{T}, \bar{\rho}),\ \ \  (\color{red}{\text{eq 2}})
         * \f}
         * The first factor \f$ \bar{\mu_0} \f$ of the product represents the viscosity in the zero-density limit and is a function of temperature only. The second factor \f$ \bar{\mu_1} \f$ represents the contribution to viscosity due to increasing density, while the third factor \f$ \bar{\mu_2} \f$ represents an enhancement of the viscosity near the critical point. The determination of each of these contributions will be consid- ered in the following sections. 
         * 
         * \f{equation}\left\{
         * \begin{matrix}
         * \bar{\mu_0}(\bar T) & =  & \frac{100\sqrt{\bar T}}{\sum\limits_{i=0}^3\frac{H_i}{\bar{T^i}}} &,& (\color{red}{\text{eq 11}}) \\
         * \bar{\mu_1}(\bar{T}, \bar{\rho})  & =  & exp\left[ \bar{\rho}\sum\limits_{i=0}^5\left(\frac{1}{\bar T} - 1\right)^i\sum\limits_{j=0}^6H_{ij}(\bar{\rho} - 1)^j \right] &,& (\color{red}{\text{eq 12}}) \\
         * \bar{\mu_2}(\bar{T}, \bar{\rho})  & = & e^{\chi_{\mu} Y}, \chi_{\mu} = 0.068 &,& (\color{red}{\text{eq 33}})
         * \end{matrix}\right.
         * \f}
         * 
         * \f{equation}
         * \xi = \xi_0\left( \frac{\bar{\chi}}{\Gamma_0}\right)^{\nu/\gamma}, \ \ \  (\color{red}{\text{eq 16}})
         * \f}
         * \f{equation}\left\{ \begin{matrix}
         * Y & = & \frac{1}{5}\color{blue}{q_C\xi}(q_D\xi)^5\left(1 - \color{blue}{q_C\xi} + (\color{blue}{q_C\xi})^2 - \frac{765}{504}(q_D\xi)^2 \right) &,& 0\le\xi\le 0.3817016416&,& (\color{red}{\text{eq 20, page 114 and figure 9}}) \\
         * Y & = & \frac{1}{12}sin(3\psi_D) - \frac{1}{4\color{blue}{q_C\xi}}sin(2\psi_D) + \frac{1}{(\color{blue}{q_C\xi})^2}\left[1- \frac{5}{4}(\color{blue}{q_C\xi})^2 \right]sin(\psi_D) - \frac{1}{(\color{blue}{q_C\xi})^3}\left\{\left[ 1 - \frac{3}{2}(\color{blue}{q_C\xi})^2 \right]\psi_D - |(\color{blue}{q_C\xi})^2-1|^{3/2}L(\omega) \right\} &,& else &,& (\color{red}{\text{eq 26}})
         * \end{matrix}\right.
         * \f}
         * \sa H2O::Table235
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param Rho Density [\f$ kg/m^3 \f$]
         * @return Dynamic viscosity [\f$ Pa s \f$]  
         */
        double mu_T_Rho(double T, double Rho);
        /**
         * @brief Dynamic viscosity
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [\f$ bar \f$]
         * @return Dynamic viscosity [\f$ Pa s \f$] 
         */
        double mu(double T, double P);
    public:
        
    };
    
}

#endif