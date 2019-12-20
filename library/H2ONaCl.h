#ifndef H2ONACL_H
#define H2ONACL_H
#include<iostream>
#include<map>  
#include<string> 
#include<cmath>
using namespace std;

// const value define
#define TMIN 273.15
#define TMAX 1273.15
#define PMIN 5e5
#define PMAX 1000e5
#define XMIN 1e-5
#define XMAX 1
// ------------
#define M_H2O   0.018015  // molar mass of water  [kg/mol]
#define M_NaCl  0.058443 // molar mass of salt(NaCl)   [kg/mol]



// define index of region
enum PhaseRegion {SinglePhase_L, TwoPhase_L_V_X0, SinglePhase_V, 
                  TwoPhase_L_H, TwoPhase_V_H, ThreePhase_V_L_H, TwoPhase_V_L_L, TwoPhase_V_L_V};
typedef map<int,string> MAP_PHASE_REGION;

struct PROP_H2ONaCl
{
    PhaseRegion Region;
    double T, H, Rho;
    double Rho_l, Rho_v, Rho_h;
    double H_l, H_v, H_h;
    double S_l, S_v, S_h;
    double X_l, X_v; // volume fraction of NaCl in vaper and liquid, it is a composition fraction. H2O + NaCl
    double Mu_l, Mu_v;
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
    double m_P, m_T, m_X; //P: Pa; T: C; X: (0, 1]
    const double *m_Parray;
    const double *m_Tarray;
    const double *m_Xarray;
    int m_number; //how man points of (P, T, X)
    PROP_H2ONaCl m_prop;
    Cr_STRUCT m_Cr;
    Cr_STRUCT init_Cr();
    f_STRUCT m_f;
    f_STRUCT init_f();
private:
    PhaseRegion findRegion(const double T, const double P, const double X, double& Xl_all, double& Xv_all);
    void fluidProp_crit_T(double T, double tol, double& P,double& Rho_l, double& Rho_v, double& h_l, double& h_v);
    void fluidProp_crit_P(double P, double tol, double& T_2ph, double& Rho_l, double& h_l, double& h_v, double& dpd_l, double& dpd_v, double& Rho_v, double& Mu_l, double& Mu_v);
public:
    cH2ONaCl(double P, double T, double X);
    cH2ONaCl(const vector<double> P, const vector<double> T, const vector<double> X);//P, T, X has same size of n
    ~cH2ONaCl();
    MAP_PHASE_REGION m_phaseRegion_name;
    void Calculate();

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
    inline double sum_array1d(double* a, int n)
    {
        double sum=0;
        for(int i=0;i<n;i++)sum+=a[i];
        return sum;
    };
    void approx_ps(double T, double& psa, double& dpsdt);
    // void two_phase_prop(double T,double Rho_l, double Rho_v, MP_STRUCT MP,ID_STRUCT ID, TWOPHASEPROP_STRUCT& l_prop,TWOPHASEPROP_STRUCT& v_prop, double& dpT);
};

#endif