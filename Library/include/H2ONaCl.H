#ifndef H2ONACL_H
#define H2ONACL_H
#include<iostream>
#include<fstream>
#include<string> 
#include<cmath>
#include<vector>
using namespace std;

// #include "omp.h"
// chose use which library to calculate water properties
// if H2ONaCl class is called by swift on ios platform, we change to use PROST to calculate water properties
// if on desktop platform, please comment this line
#include "stdfunc.H"
// #include "Fluid.H"
#include "H2O.H"
#include "NaCl.H"
#include "dataStruct_H2ONaCl.H"

#define USE_PROST 1  

#ifdef USE_PROST
    #include "steam4.h"
#else 
	#include "IAPWS-97.H"
#endif 

namespace H2ONaCl
{
    // ============= Constants of H2O-NaCl =======================================
    double const PMIN = 1; /**< Minimum valid pressure of H2O, [bar]. */ 
    double const PMAX = 1000; /**< Minimum valid pressure of H2O, [bar]. */ 
    double const TMIN = 273.15; /**< Minimum valid pressure of H2O, [K].  */ //TODO： 最后统一为C
    double const TMAX = 1273.15; /**< Minimum valid pressure of H2O, [K]. */ //TODO: 最后统一为C
    double const TMIN_K = 273.15; /**< Minimum valid pressure of H2O, [K].  */ 
    double const TMAX_K = 1273.15; /**< Minimum valid pressure of H2O, [K]. */ 
    double const TMIN_C = 0; /**< Minimum valid pressure of H2O, [C].  */ 
    double const TMAX_C = 1000; /**< Minimum valid pressure of H2O, [C]. */ 
    double const XMIN = 0;
    double const XMAX = 1;
    // ===========================================================================

    /**
     * @brief EOS and thermodynamic properties of \f$H_2O-NaCl\f$ system. 
     * 
     * This implementation is based published references \cite Driesner2007Part1, \cite Driesner2007Part2 and \cite Klyukin2020. 
     * 
     */
    class cH2ONaCl
    {
    private:
        void init_PhaseRegionName();
        void init_prop();
        double m_P, m_T, m_Xwt, m_Xmol; //P: Pa; T: C; X, wt%: (0, 1]. 
        // Note that !!! m_T unite is C, but to keep consistent with OpenFoam, the T variable in public member function with unit of K
        const double *m_Parray;
        const double *m_Tarray;
        const double *m_Xarray;
        int m_number; //how many points of (P, T, X)
        Cr_STRUCT m_Cr;
        Cr_STRUCT init_Cr();
        f_STRUCT m_f;
        f_STRUCT init_f();
        bool m_colorPrint;
    private:
        PhaseRegion findRegion(const double T, const double P, const double X, double& Xl_all, double& Xv_all);
        void calcRho(int reg, double T_in, double P_in, double X_l, double X_v, double& Rho_l, double& Rho_v, double& Rho_h, 
                        double& V_l_out, double& V_v_out, double& T_star_l_out, double& T_star_v_out, double& n1_v_out, double& n2_v_out);
        void calcEnthalpy(int reg, double T_in, double P_in, double X_l, double X_v,
            double& h_l, double& h_v, double& h_h);
        void calcViscosity(int reg, double P, double T, double Xw_l, double Xw_v, double& mu_l, double& mu_v);
        void fluidProp_crit_T(double T, double tol, double& P,double& Rho_l, double& Rho_v, double& h_l, double& h_v);
        void fluidProp_crit_P(double P, double tol, double& T_2ph, double& Rho_l, double& h_l, double& h_v, double& dpd_l, double& dpd_v, double& Rho_v, double& Mu_l, double& Mu_v);
        void guess_T_PhX(double P, double h, double X, double& T1, double& T2);
        void calc_sat_lvh(PROP_H2ONaCl& prop, double h, double X, bool isDeriv=false);
        void calc_halit_liqidus(double Pres, double Temp, double& X_hal_liq, double& T_hm);
        // function of water properties: using freesteam or PROST, if both of them are available, default to use freesteam
        double water_rho_pT(double p, double T_K);
        double water_h_pT(double p, double T_K);
        double water_mu_pT(double p, double T_K);
    public:
        // cH2ONaCl(double P, double T_K, double X);//P: Pa. T: K  X, wt%: (0, 1]
        cH2ONaCl();
        // cH2ONaCl(const std::vector<double> P, const std::vector<double> T, const std::vector<double> X);//P, T, X has same size of n
        ~cH2ONaCl();
        H2O::cH2O m_water;
        NaCl::cNaCl m_NaCl;
        MAP_PHASE_REGION m_phaseRegion_name;
        void prop_pTX(double p, double T_K, double X_wt, bool visc_on=true);
        void prop_pHX(double p, double H, double X_wt); /** Calculate properties by P, H, X */
        double rho_pTX(double p, double T_K, double X_wt); //get bulk density. p: Pa; T: K; X: wt%
        double rho_l_pTX(double p, double T_K, double X_wt); //get density of liquid. p: Pa; T: K; X: wt%
        double mu_l_pTX(double p, double T_K, double X_wt); //get dynamic viscosity of liquid. p: Pa; T: K; X: wt%
        double mu_pTX(double p, double T_K, double X_wt); //get bulk dynamic viscosity. p: Pa; T: K; X: wt%
        PROP_H2ONaCl m_prop;
        //X:[0,1]; T: deg. C; P: bar
        void writeProps2VTK(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<H2ONaCl::PROP_H2ONaCl> props, std::string fname, bool isWritePy=true, std::string xTitle="x", std::string yTitle="y", std::string zTitle="z");
        void writeProps2xyz(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<H2ONaCl::PROP_H2ONaCl> props, std::string fname, std::string xTitle="x", std::string yTitle="y", std::string zTitle="z", string delimiter=" ");
        friend ostream & operator<<(ostream & out,  cH2ONaCl & A);
        void setColorPrint(bool colorPrint){m_colorPrint=colorPrint;}
        string checkTemperatureRange(double temperature_C);
        string checkPressureRange(double pressure_bar);
        string checkSalinityRange(double salinity);
        string CheckRange_H(double H0, double P0, double X0);
        string CheckRanges_H_PX(double HMIN0, double HMAX0, double PXrange[4]);
        string CheckRanges_H_P(double HMIN0, double HMAX0, double Prange[2], double X0);
        /**
         * @brief Check valid ranges of enthalpy and salinity
         * 
         * @param HMIN0 
         * @param HMAX0 
         * @param Xrange 
         * @param P0 
         * @return string 
         */
        string CheckRanges_H_X(double HMIN0, double HMAX0, double Xrange[2], double P0);

        /**
         * @brief Critical curve of NaCl-H2O system. It can be evaluated into three segments, see equation (5) - (7) of reference \cite Driesner2007Part1. 
         * 
         * 
         * @param T Input temperature, [\f$ ^{\circ}\text{C} \f$]
         * @param P_crit Output critical pressure, [bar]
         * @param X_crit Output critical salinity, [Mole fraction of NaCl]
         * 
         */
        void P_X_Critical(double T, double& P_crit, double& X_crit);
        /**
         * @brief The composition of halite-saturated liquid, \f$ X_{NaCl, sat}^L \f$ (the halite liquidus), is rather well known at high temperatures from about 400 \f$^{\circ}C \f$ to the melting curve of NaCl at pressures to about 4000 bar. See equation (8) of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_HaliteLiquidus(double T, double P);
        /**
         * @brief Composition of halite-saturated vapor. See equation (9) and (14) of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_VaporHaliteCoexist(double T, double P);
        /**
         * @brief Temperature-pressure relations when vapor-liquid-halite coexist. See equation (10) and Table 6 of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @return double Pressure [bar]
         */
        double P_VaporLiquidHaliteCoexist(double T);
        /**
         * @brief Salinity on liquid branch of Vapor + Liquid coexist surface. See equation (11) and Table 7 of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_VaporLiquidCoexistSurface_LiquidBranch(double T, double P);
        /**
         * @brief Salinity on vapor branch of Vapor + Liquid coexist surface. See equation (13-17) and Table 8 of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_VaporLiquidCoexistSurface_VaporBranch(double T, double P);
        /**
         * @brief \f$ T_V^* \f$ is used to determin the correction between NaCl-H2O and H2O. For a given molar volume, NaCl solution and pure water temperature are \f$ T \f$ and \f$ T_V^* \f$, respectively. The principle is shown in Fig.2 (below) of reference \cite Driesner2007Part2. \f$ T_V^* \f$ is calculated from equation (8, 13) of \cite Driesner2007Part2.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @param X Salinity [Mole fraction of NaCl]
         * @return double Temperature [\f$ ^{\circ}\text{C} \f$]
         */
        double T_star_V(double T, double P, double X);
        /**
         * @brief Extrapolated molar volume outside range of correlation. See equation (17) of reference \cite Driesner2007Part2.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @param X Salinity [Mole fraction of NaCl]
         * @return double \f$ V_{extrapol} \f$ [\f$cm^3\ mol^{-1} \f$] 
         */
        double V_extrapol(double T, double P, double X);
        inline double Wt2Mol(double X_wt)
        {
            return X_wt/NaCl::MolarMass/(X_wt/NaCl::MolarMass + (1-X_wt)/H2O::MolarMass);
        };
        /**
         * @brief Only used in #V_extrapol
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @param X Salinity [Mole fraction of NaCl]
         * @return double Density [\f$ kg/m^3 \f$]
         */
        double Rho_Br_for_V_extrapol(double T, double P, double X);
        /**
         * @brief Density of brine.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @param X Salinity [Mole fraction of NaCl]
         * @return double Density [\f$ kg/m^3 \f$]
         */
        double Rho(double T, double P, double X);
    private:
        
        inline double Xwt2Xmol(double X){return (X/NaCl::MolarMass)/(X/NaCl::MolarMass+(1-X)/H2O::MolarMass);};
        void approx_Rho_lv(double T, double& Rho_l , double& Rho_v);
        MP_STRUCT bb(double T);
        ID_STRUCT ideal(double T);
        void twoPhaseProp(double T,double Rho_l, double Rho_v, MP_STRUCT MP,ID_STRUCT ID, TWOPHASEPROP_STRUCT& l_prop, TWOPHASEPROP_STRUCT& v_prop);
        void psatc(double T, double& h_l, double& h_v ,double& Rho_l, double& Rho_v,double& dpd_l,double& dpd_v, double& psa);
        BS_STRUCT base(double T, double Rho_l, MP_STRUCT MP);
        RS_STRUCT resid(double T, double Rho_l);
        TWOPHASEPROP_STRUCT props(double T, double Rho, BS_STRUCT BS, RS_STRUCT RS, ID_STRUCT ID);
        void approx_ps(double T, double& psa, double& dpsdt);
        // double water_tp_IAPS84(double P, double T, double Rho, double& dRhodP, double& h, double& Mu, double tol=1e9,bool validation=true);
        // int region_tp(double T, double P);
        template <typename T> T sum_array1d(T* a, int n);
        template <typename T> T max(std::vector<T> data);
        template <typename T> T min(std::vector<T> data);
    };
}

#endif