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
enum PhaseRegion {SinglePhase_L, SinglePhase_V, SinglePhase_H, TwoPhase_V_L,TwoPhase_V_H, TwoPhase_L_H, ThreePhase_V_L_H};
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
private:
    PhaseRegion findRegion(const double T, const double P, const double X);
    void fluidProp_crit_T(double T, double tol, double& P,double& Rho_l, double& Rho_v, double& h_l, double& h_v);
public:
    cH2ONaCl(double P, double T, double X);
    // cH2ONaCl(const double* P, const double* T, const double* X, const int n);//P, T, X has same size of n
    ~cH2ONaCl();
    MAP_PHASE_REGION m_phaseRegion_name;
    void Calculate();

private:
    inline double Xwt2Xmol(double X){return (X/M_NaCl)/(X/M_NaCl+(1-X)/M_H2O);};
    void approx_Rho_lv(double T, double& Rho_l , double& Rho_v);
};

#endif