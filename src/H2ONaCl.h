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

//  #define PLATFORM_IOS 1  

#ifdef PLATFORM_IOS
    #include "steam4.h"
#else 
	#include "IAPWS-97.H"
#endif 
#include "stdio.h"

#include<map>  
namespace SWEOS
{
    #ifdef _WIN32
        // define color, this seems only work on MacOS and linux, doesn't work on windows
        #define ERROR_COUT ""
        #define WARN_COUT ""
        #define COLOR_PURPLE ""
        #define COLOR_RED ""
        #define COLOR_GREEN ""
        #define COLOR_YELLOW ""
        #define COLOR_BLUE ""
        #define COLOR_DEFAULT ""
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
    // const value define
    double const TMIN=273.15;
    double const TMAX=1273.15;
    double const PMIN=1e5;
    double const PMAX = 1000e5;
    double const XMIN = 0;
    double const XMAX = 1;
    // --------------
    double const Kelvin= 273.15;
    // ------------
    double const M_H2O = 0.018015; // molar mass of water  [kg/mol]
    double const M_NaCl = 0.058443; // molar mass of salt(NaCl)   [kg/mol]

    // double const P_crit =2.2054915e7;


    // define index of region
    enum PhaseRegion {SinglePhase_L, TwoPhase_L_V_X0, SinglePhase_V, 
                    TwoPhase_L_H, TwoPhase_V_H, ThreePhase_V_L_H, TwoPhase_V_L_L, TwoPhase_V_L_V};
                    
    typedef std::map<int,std::string> MAP_PHASE_REGION;

    struct PROP_H2ONaCl
    {
        PhaseRegion Region;
        double T, H, Rho, Mu; //temperature, bulk enthalpy, bulk density
        double Rho_l, Rho_v, Rho_h; //density
        double H_l, H_v, H_h; //enthalpy
        double S_l, S_v, S_h; //saturation
        double X_l, X_v; // volume fraction of NaCl in vaper and liquid, it is a composition fraction. H2O + NaCl
        double Mu_l, Mu_v;//viscosity
        // derivertive
        // double dS_hdh, dS_vdh, dS_ldh, dRhodh;
    };

    struct MP_STRUCT
    {
        double b1, b1t, b1tt;
        double b2, b2t, b2tt;
    };
    struct ID_STRUCT
    {
        double f, ft, ftt;
    };
    struct TWOPHASEPROP_STRUCT
    {
        double f, p, s, g, u, h, dpd, dpt, cv, x;
    };
    struct BS_STRUCT
    {
        double x, f;
        double fd, fdd, ft, ftd, ftt;
    };
    struct RS_STRUCT
    {
        double f, ft, ftd;
        double fd, fdd, ftt;
    };
    struct Cr_STRUCT
    {
        double g[9][6], k[4], l[4], gg[4], t[4], d[4], a[4], b[4];
    };
    struct f_STRUCT
    {
        double f[11];
        double sum_f10;
    };


    /**
     * @brief Class of thermal model of seawater
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
        MAP_PHASE_REGION m_phaseRegion_name;
        void prop_pTX(double p, double T_K, double X_wt, bool visc_on=true);
        void prop_pHX(double p, double H, double X_wt); /** Calculate properties by P, H, X */
        double rho_pTX(double p, double T_K, double X_wt); //get bulk density. p: Pa; T: K; X: wt%
        double rho_l_pTX(double p, double T_K, double X_wt); //get density of liquid. p: Pa; T: K; X: wt%
        double mu_l_pTX(double p, double T_K, double X_wt); //get dynamic viscosity of liquid. p: Pa; T: K; X: wt%
        double mu_pTX(double p, double T_K, double X_wt); //get bulk dynamic viscosity. p: Pa; T: K; X: wt%
        PROP_H2ONaCl m_prop;
        //X:[0,1]; T: deg. C; P: bar
        void writeProps2VTK(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<PROP_H2ONaCl> props, std::string fname, bool isWritePy=true, string xTitle="x", string yTitle="y", string zTitle="z");
        friend ostream & operator<<(ostream & out,  cH2ONaCl & A);
        void setColorPrint(bool colorPrint){m_colorPrint=colorPrint;}
        string checkTemperatureRange(double temperature_C);
        string checkPressureRange(double pressure_bar);
        string checkSalinityRange(double salinity);
        string CheckRange_H(double H0, double P0, double X0);
        string CheckRanges_H_PX(double HMIN0, double HMAX0, double PXrange[4]);
        string CheckRanges_H_P(double HMIN0, double HMAX0, double Prange[2], double X0);
        string CheckRanges_H_X(double HMIN0, double HMAX0, double Xrange[2], double P0);
    private:
        inline double Xwt2Xmol(double X){return (X/M_NaCl)/(X/M_NaCl+(1-X)/M_H2O);};
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