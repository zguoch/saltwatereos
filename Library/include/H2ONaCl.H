#ifndef H2ONACL_H
#define H2ONACL_H
#include<iostream>
#include<fstream>
#include<string> 
#include<cmath>
#include<vector>
using namespace std;
#include <stdlib.h>
// #include "omp.h"
// chose use which library to calculate water properties
// if H2ONaCl class is called by swift on ios platform, we change to use PROST to calculate water properties
// if on desktop platform, please comment this line
#include "stdfunc.H"
// solve roots of general polynomial
// This code comes from https://www.codeproject.com/Articles/674149/A-Real-Polynomial-Class-with-Root-Finder
#include "Polynomial.h"

// #include "Fluid.H"
#include "H2O.H"
#include "NaCl.H"
#include "dataStruct_H2ONaCl.H"
#include "LookUpTableForest.h"

#if USE_OMP == 1
#include "omp.h"
#endif

#define USE_PROST 1  

#ifdef USE_PROST
    #include "steam4.h"
#else 
	#include "IAPWS-97.H"
#endif 

namespace H2ONaCl
{
    // ============= Constants of H2O-NaCl =======================================
    double const PMIN = 5; /**< Minimum valid pressure of H2O-NaCl, [bar]. */ 
    double const PMAX = 5000; /**< Maximum valid pressure of H2O-NaCl, [bar]. */ 
    double const TMIN = 273.25; /**< Minimum valid pressure of H2O-NaCl, [K].  */ 
    double const TMAX = 1273.15; /**< Maximum valid pressure of H2O-NaCl, [K]. */
    double const TMIN_K = 273.25; /**< Minimum valid pressure of H2O-NaCl, [K].  */ 
    double const TMAX_K = 1273.15; /**< Maximum valid pressure of H2O-NaCl, [K]. */ 
    double const TMIN_C = 0.1; /**< Minimum valid temperature of H2O-NaCl, [C].  */ 
    double const TMAX_C = 1000; /**< Maximum valid pressure of H2O-NaCl, [C]. */ 
    double const XMIN = 0;  /**< Minimum valid salinity, [mol fraction].  */ 
    double const XMAX = 1;  /**< Maximum valid salinity, [mol fraction].  */ 
    double const HMAX = 6E6; /**< Approximately maximum valid enthalpy, [J/kg]. It may be used to normize vtk output file data*/
    // ===========================================================================
    double const P_max_LVH  = 390.14744433796; /**< Maximum pressure of L+V+H region, [bar] */
    double const T_Pmax_LVH = 594.63244000000; /**< Temperature at that P_max_LVH */

    // ======
    template <int dim>
    struct FIELD_DATA
    {
        // point data field
        // H2ONaCl::PROP_H2ONaCl prop_point[1<<dim]; // properties at vertiex
        // cell data field
        LOOKUPTABLE_FOREST::NeedRefine need_refine; // indicator of what kind of the need-refined quad position (phase boundary), or what kind of properties need to refine
        H2ONaCl::PhaseRegion phaseRegion_cell;
    };
    // ============= bitmask of properties =======================================
    #define Update_prop_rho     2   // 2^1 = 2
    #define Update_prop_h       4   // 2^2 = 4
    #define Update_prop_drhodh  8   // 2^3 = 8
    #define Update_prop_T       16   // 2^4 = 16
    // ===========================================================================
    typedef LOOKUPTABLE_FOREST::LookUpTableForest<2, H2ONaCl::FIELD_DATA<2> > LookUpTableForest_2D;
    typedef LOOKUPTABLE_FOREST::LookUpTableForest<3, H2ONaCl::FIELD_DATA<3> > LookUpTableForest_3D;

    /**
     * @brief Check the criterion of a property and determine if need to refine. Should use relative error criterion, rather absolute error otherwise it is not fair for vapor region.
     * 
     */
    #define CHECK_REFINE_PROP_RMSD(PROP) \
    { \
        double mean_Prop = props_refine_check[0].PROP; \
        for(int i=0;i<forest->m_num_children;i++)mean_Prop += props_refine_check[i+1].PROP; \
        mean_Prop = mean_Prop / (forest->m_num_children + 1);  \
        double RMSD_Prop = pow((props_refine_check[0].PROP - mean_Prop)/mean_Prop, 2.0); \
        for(int i=0;i<forest->m_num_children;i++)RMSD_Prop += pow((props_refine_check[i+1].PROP - mean_Prop)/mean_Prop, 2.0); \
        RMSD_Prop = sqrt(RMSD_Prop/(forest->m_num_children + 1)); \
        if(RMSD_Prop > forest->m_RMSD_RefineCriterion.PROP) \
        { \
            data->need_refine = LOOKUPTABLE_FOREST::NeedRefine_##PROP; \
        } \
    }

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
        void init_prop(H2ONaCl::PROP_H2ONaCl& prop);
        Cr_STRUCT m_Cr;
        Cr_STRUCT init_Cr();
        f_STRUCT m_f;
        f_STRUCT init_f();
        bool m_colorPrint;
        TABLE4 m_tab4_Driesner2007a; /**< Parameters for the critical curve Critical */
    private:
        void createTable4_Driesner2007a(TABLE4& table4);
        /**
         * @brief Find phase region when (T,P,X) is given.
         * 
         * @param T Temperature [deg.C]
         * @param P Pressure [Pa]
         * @param X Bulk salinity [Mole fraction of NaCl]
         * @param Xl_all Salinity of liquid phase [mole fraction of NaCl]
         * @param Xv_all Salinity of vapor phase [mole fraction of NaCl]
         * @return H2ONaCl::PhaseRegion 
         */
        PhaseRegion findRegion(const double T, const double P, const double X, double& Xl_all, double& Xv_all);
        
        void calcRho(int reg, double T_in, double P_in, double X_l, double X_v, double& Rho_l, double& Rho_v, double& Rho_h, 
                        double& V_l_out, double& V_v_out, double& T_star_l_out, double& T_star_v_out, double& n1_v_out, double& n2_v_out);
        /**
         * @brief Calculate specific enthalpy.
         * 
         * @param reg Phase region index, H2ONaCl::PhaseRegion
         * @param T_in Temperature [deg.C]
         * @param P_in Pressure [bar]
         * @param X_l Salinity of liquid phase [mole fraction of NaCl]
         * @param X_v Salinity of vapor phase [mole fraction of NaCl]
         * @param h_l Specific enthalpy of liquid phase [J/kg]
         * @param h_v Specific enthalpy of vapor phase [J/kg]
         * @param h_h Specific enthalpy of halite [J/kg]
         */
        void calcEnthalpy(int reg, double T_in, double P_in, double X_l, double X_v,
            double& h_l, double& h_v, double& h_h);
        /**
         * @brief Calculate dynamic viscosity.
         * 
         * @param reg Phase region index, H2ONaCl::PhaseRegion
         * @param P Pressure [Pa]
         * @param T Temperature [deg.C]
         * @param Xw_l \todo 搞清楚这里的盐度单位和黏度计算公式
         * @param Xw_v 
         * @param mu_l 
         * @param mu_v 
         */
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
        H2ONaCl::MAP_PHASE_REGION m_phaseRegion_name;
        H2ONaCl::PROP_H2ONaCl m_prop;
    // ========================== benchmark test have been done ====================
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
         * @brief Array version.
         * 
         * @param T 
         * @param P_crit 
         * @param X_crit 
         */
        void P_X_Critical(std::vector<double> T, std::vector<double>& P_crit, std::vector<double>& X_crit);
    // ==================================================================================================================

        inline string getPhaseRegionName(PhaseRegion regionID){return m_phaseRegion_name[regionID];};
        /**
         * @brief Calculate thermal dynamic properties of NaCl-H2O system.
         * 
         * @param p pressure [Pa]
         * @param T_K Temperature [K]
         * @param X_wt Salinity [mass fraction, [0,1]]
         * @param visc_on Calculate viscosity or not, default is true.
         * @return H2ONaCl::PROP_H2ONaCl 
         */
        H2ONaCl::PROP_H2ONaCl prop_pTX(double p, double T_K, double X_wt, bool visc_on=true);
       
        /**
         * @brief Calculate thermal dynamic properties of NaCl-H2O system.
         * 
         * @param p pressure [Pa]
         * @param H specific enthalpy [J/kg]
         * @param X_wt Salinity [mass fraction, [0,1]]
         * @return H2ONaCl::PROP_H2ONaCl 
         */
        H2ONaCl::PROP_H2ONaCl prop_pHX(double p, double H, double X_wt);
        /**
         * @brief Calculate bulk density.
         * 
         * @param p [Pa]
         * @param T_K [K]
         * @param X_wt [mass fraction of NaCl]
         * @return double Density [kg/m3]
         */
        double rho_pTX(double p_Pa, double T_K, double X_wt); 
        double rho_l_pTX(double p_Pa, double T_K, double X_wt); 
        double mu_l_pTX(double p_Pa, double T_K, double X_wt);
        double mu_pTX(double p_Pa, double T_K, double X_wt);
        void writeProps2VTK(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<H2ONaCl::PROP_H2ONaCl> props, std::string fname, bool isNormalize=true, std::string xTitle="x", std::string yTitle="y", std::string zTitle="z");
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
         * @brief Check valid ranges of enthalpy given salinity
         * 
         * @param HMIN0 
         * @param HMAX0 
         * @param Xrange 
         * @param P0 
         * @return string 
         */
        string CheckRanges_H_X(double HMIN0, double HMAX0, double Xrange[2], double P0);

        /**
         * @brief Calculate critical T given P
         * 
         * @param P Input pressure, [bar]
         * @param T_crit Output critical temperature, [\f$ ^{\circ}\text{C} \f$]
         * @param X_crit Output critical salinity, [Mole fraction of NaCl]
         * @return double 
         */
        void T_X_Critical(double P, double& T_crit, double& X_crit);
        /**
         * @brief Calculate critical curve of \f$H_2O-NaCl \f$ system in the whole valid region, and then write as to file in format of VTK or dat.
         * 
         * @param filename File name without extension, default filename is "CriticalCurve", i.e. the default output path is the current path.
         * @param Tmin Minimum temperature (in unit of \f$ ^{\circ}C\f$) of temperature for critical curve, default is TMIN_C.
         * @param Tmax Maximum temperature (in unit of \f$ ^{\circ}C\f$) of temperature for critical curve, default is TMIN_C.
         * @param dT Temperature spacing of the critical curve. Default is 1 \f$ ^{\circ}C\f$.
         * @param fmt Output file format. Option is one of H2ONaCl::fmt_vtk, H2ONaCl::fmt_dat.
         */
        void writeCriticalCurve(string filename="CriticalCurve",double Tmin=H2O::T_Critic, double Tmax=TMAX_C, double dT=1, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk);
        void writeNaClMeltingCurve(string filename="NaClMeltingCurve",double Pmin=PMIN, double Pmax=PMAX, double dP=1, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk);
        void writeVaporLiquidHalite_V_L_H_Curve(string filename="VaporLiquidHalite",double Tmin=TMIN_C, double Tmax=NaCl::T_Triple, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk, int nT=100);
        void writeH2OBoilingCurve(string filename="BoilingCurve_H2O",double Tmin=TMIN_C, double Tmax=H2O::T_Critic, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk, int nT=100);
        /**
         * @brief The composition of halite-saturated liquid, \f$ X_{NaCl, sat}^L \f$ (the halite liquidus), is rather well known at high temperatures from about 400 \f$^{\circ}C \f$ to the melting curve of NaCl at pressures to about 4000 bar. See equation (8) of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_HaliteLiquidus(double T, double P);
        std::vector<double> X_HaliteLiquidus(std::vector<double> T, std::vector<double> P);
        /**
         * @brief Write halite liquidus surface as vtk file format
         * 
         * @param filename 
         * @param Tmin [C]
         * @param Tmax [C]

         * @param Pmax [bar]
         * @param fmt 
         * @param nT 
         * @param nP 
         */
        void writeHaliteLiquidusSurface(string filename="HaliteLiquidus",double Tmin=TMIN_C, double Tmax=NaCl::T_Triple, double Pmax=PMAX, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk, int nT=100, int nP=200);
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
         * @brief Get halite and salinity of V+L+H boundary surface.
         * 
         * Of course this is only valid for pressure in range of [H2ONaCl::PMIN, H2ONaCl::cH2ONaCl::Pmax_VaporLiquidHaliteCoexist (~390.147 bar)]
         * 
         * @param P Pressure [bar]
         * @return  std::vector<double> [Hmin, Hmax, Xmin(wt), Xmax(wt)] with four components, or empty for invalid P
         */
        std::vector<double> HX_VaporLiquidHaliteCoexist(double P);
        /**
         * @brief Get maximum pressure and corresponding temperature of V+L+H surface
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return  
         */
        void Pmax_VaporLiquidHaliteCoexist(double& T, double& P);
        /**
         * @brief Solve temperature of vapor-liquid-halite coexist boundary when given P
         * 
         * @param P Pressure [bar]
         * @return double Temperature [\f$ ^{\circ}\text{C} \f$]
         */
        std::vector<double> T_VaporLiquidHaliteCoexist(double P);
        /**
         * @brief Write vapor+Liquid+Halite coexist surface
         * 
         * @param Tmin [C]
         * @param Tmax [C]
         * @param dT [C]
         */
        void writeVaporLiquidHaliteCoexistSurface(string filename="VaporLiquidHalite", double Tmin=TMIN_C, double Tmax=NaCl::T_Triple, double dT=1, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk);
        /**
         * @brief Salinity on liquid branch of Vapor + Liquid coexist surface. See equation (11) and Table 7 of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_VaporLiquidCoexistSurface_LiquidBranch(double T, double P);
        std::vector<double> X_VaporLiquidCoexistSurface_LiquidBranch(std::vector<double> T, std::vector<double> P);
        /**
         * @brief Salinity on vapor branch of Vapor + Liquid coexist surface. See equation (13-17) and Table 8 of reference \cite Driesner2007Part1.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @return double Salinity [Mole fraction of NaCl]
         */
        double X_VaporLiquidCoexistSurface_VaporBranch(double T, double P);
        std::vector<double> X_VaporLiquidCoexistSurface_VaporBranch(std::vector<double> T, std::vector<double> P);

        void writeVaporLiquidCoexistSurface(string filename="VaporLiquid", double Tmin=TMIN_C, double Tmax=TMAX_C, H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk, int nT=100, int nP=100); 
        void T_star_V_n1n2(double P, double X, double& n1, double& n2);
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

        /**
         * @brief Convert mass fraction of NaCl to molar fraction. 
         * 
         * @param X_wt [0,1]
         * @return double [0,1]
         */
        inline double Wt2Mol(double X_wt)
        {
            return X_wt/NaCl::MolarMass/(X_wt/NaCl::MolarMass + (1-X_wt)/H2O::MolarMass);
        };
        /**
         * @brief Convert molar fraction of NaCl to mass fraction.
         * 
         * @param X_mol [0,1]
         * @return double [0,1]
         */
        inline double Mol2Wt(double X_mol)
        {
            return NaCl::MolarMass * X_mol / (NaCl::MolarMass * X_mol + (1 - X_mol) * H2O::MolarMass);
        };
        inline std::vector<double> Mol2Wt(std::vector<double> X_mol){
            std::vector<double> x_wt;
            for (size_t i = 0; i < X_mol.size(); i++)
            {
                x_wt.push_back(Mol2Wt(X_mol[i]));
            }
            return x_wt;
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
        double Rho_brine(double T, double P, double X);

        /**
         * @brief Find phase region for giving T, P and X.
         * 
         * @param T Temperature [\f$ ^{\circ}\text{C} \f$]
         * @param P Pressure [bar]
         * @param X_wt Salinity [Mass fraction of NaCl]
         * @param Xl_all Salinity of liquid phase [mole fraction of NaCl]
         * @param Xv_all Salinity of vapor phase [mole fraction of NaCl]
         * @return H2ONaCl::PhaseRegion 
         */
        PhaseRegion findPhaseRegion(const double T_c, const double P_bar, const double X_wt, double& Xl_all, double& Xv_all);
        PhaseRegion findPhaseRegion(const double T_c, const double P_bar, const double X_wt);
        /**
         * @brief This function has the same variable order and unit as prop_pTX.
         * 
         * @param p_Pa
         * @param T_K 
         * @param X_wt 
         * @return PhaseRegion 
         */
        PhaseRegion findPhaseRegion_pTX(double p_Pa, double T_K, double X_wt);
        /**
         * @brief Write VLH phase boundary surface in PHX space to VTU file
         * 
         * @param scale_X scale of X, default is 1
         * @param scale_H scale of H, default is 1.0/H2ONaCl::cH2ONaCl::HMAX
         * @param scale_P scale of P, default is 1.0/H2ONaCl::cH2ONaCl::PMAX
         * @param filename output vtu filename, default is "VLH.vtu"
         * @param fmt file format of the output file, default is vtk format
         * @param nP data points in P axis
         */
        void writePhaseSurface_XHP(double scale_X=1, double scale_H=1.0/HMAX, double scale_P=1.0/1000, string outpath="./", H2ONaCl::fmtOutPutFile fmt=H2ONaCl::fmt_vtk, int nP=500);
    private:
        template<int dim>
        void interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf, double* xyz_min_target, H2ONaCl::PROP_H2ONaCl& prop, const double xyz[dim]);
        
        template<int dim>
        void interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim,H2ONaCl::FIELD_DATA<dim> > *targetLeaf, double* xyz_min_target, double* props, const double xyz[dim]);
        
        void init_supported_props();
    public:
        int m_num_threads;
    #if USE_OMP == 1
        void set_num_threads(int num_threads);
        int get_num_threads();
    #endif
        int m_dim_lut;
        std::map<int, propInfo> m_supported_props;
        std::map<int, propInfo> m_update_which_props;
        void *m_pLUT;
        void parse_update_which_props(int update_which_props);
        /**
         * @brief Create a LUT 2D object in PTX space. Create different 2D LUT according to type and  xy limits, then access through member variable m_lut_PTX_2D
         * 
         * @param type 
         * @param xy_min 
         * @param xy_max 
         * @param z 
         */
        void createLUT_2D(double xy_min[2], double xy_max[2], double constZ, LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var, LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level = 4, int max_level = 6, int update_which_props=0);
        void createLUT_2D(double xmin, double xmax, double ymin, double ymax, double constZ, LOOKUPTABLE_FOREST::CONST_WHICH_VAR const_which_var, LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level = 4, int max_level = 6, int update_which_props=0);
        void createLUT_3D(double xyz_min[3], double xyz_max[3], LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level = 4, int max_level = 6, int update_which_props=0);
        void createLUT_3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, LOOKUPTABLE_FOREST::EOS_ENERGY TorH, int min_level = 4, int max_level = 6, int update_which_props=0);
        LOOKUPTABLE_FOREST::Quadrant<2,H2ONaCl::FIELD_DATA<2> > *lookup(H2ONaCl::PROP_H2ONaCl& prop, double x, double y);
        LOOKUPTABLE_FOREST::Quadrant<2,H2ONaCl::FIELD_DATA<2> > *lookup(double* props, double* xyz_min_target, double x, double y, bool is_cal=true);
        LOOKUPTABLE_FOREST::Quadrant<3,H2ONaCl::FIELD_DATA<3> > *lookup(H2ONaCl::PROP_H2ONaCl& prop, double x, double y, double z);
        LOOKUPTABLE_FOREST::Quadrant<3,H2ONaCl::FIELD_DATA<3> > *lookup(double* props, double* xyz_min_target, double x, double y, double z, bool is_cal=true);
        H2ONaCl::PROP_H2ONaCl lookup(double x, double y); //for python API
        H2ONaCl::PROP_H2ONaCl lookup(double x, double y, double z); //for python API
        LOOKUPTABLE_FOREST::Quadrant<2,H2ONaCl::FIELD_DATA<2> > *lookup_only(H2ONaCl::PROP_H2ONaCl& prop, double x, double y);
        LOOKUPTABLE_FOREST::Quadrant<3,H2ONaCl::FIELD_DATA<3> > *lookup_only(H2ONaCl::PROP_H2ONaCl& prop, double x, double y, double z);
        H2ONaCl::PROP_H2ONaCl lookup_only(double x, double y); //for python API
        H2ONaCl::PROP_H2ONaCl lookup_only(double x, double y, double z); //for python API
        void destroyLUT();
        void loadLUT(string filename);
        LookUpTableForest_2D* getLUT_2D(); //for Python API
        LookUpTableForest_3D* getLUT_3D(); //for Python API
        void save_lut_to_vtk(string filename);
        void save_lut_to_binary(string filename);
        // ========
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
        template <typename T> T min(T a, T b){return (a>b ? b : a);};
        template <typename T> T max(T a, T b){return (a<b ? b : a);};
    private:
        /**
         * @brief Write array of X,Y,Z to VTU (Serial vtkUnstructuredGrid) file as a 3D line.
         * 
         * @param filename 
         * @param X 
         * @param Y 
         * @param Z 
         * @return void 
         */
        void writeVTK_PolyLine(string filename,vector<double> X, vector<double> Y, vector<double> Z);
        void writeVTK_Triangle_Strip(string filename, vector<vector<double> > X, vector<vector<double> > Y, vector<vector<double> > Z, double scale_X=1.0, double scale_Y=1.0, double scale_Z=1.0);
        /**
         * @brief (X[i,:], Y[i,:], Z[i,:]) describe a polygon plane. ([X[i:i+1,j], X[i:i+1, j+1]]) construct a VTK_QUAD in the normal direction of polygon plane.
         * 
         * @param filename 
         * @param X N X M 2D array, M is the node number of each polygon, N is layers number.
         * @param Y N X M 2D array
         * @param Z N X M 2D array
         * @param scale_X Scale factor in X direction, default is 1.0
         * @param scale_Y Scale factor in Y direction, default is 1.0
         * @param scale_Z Scale factor in Z direction, default is 1.0
         * @param includeTwoEndsPolygon whether include polygon on two ends, default is true
         */
        void writeVTK_Quads(string filename, vector<vector<double> > X, vector<vector<double> > Y, vector<vector<double> > Z, double scale_X=1.0, double scale_Y=1.0, double scale_Z=1.0, bool includeTwoEndsPolygon=true);
    };
}

#endif