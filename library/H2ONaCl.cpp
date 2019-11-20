#include "H2ONaCl.h"

cH2ONaCl::cH2ONaCl(double P, double T, double X)
:m_P(P),
m_T(T),
m_X(Xwt2Xmol(X)),
m_Parray(NULL),
m_Tarray(NULL),
m_Xarray(NULL),
m_number(1)
{
    init_PhaseRegionName();
}
// cH2ONaCl::cH2ONaCl(const double* P, const double* T, const double* X, const int n)
// :m_P(0),
// m_T(0),
// m_X(0),
// m_Parray(P),
// m_Tarray(T),
// m_Xarray(X),
// m_number(n)
// {
//     init_PhaseRegionName();
// }

cH2ONaCl::~cH2ONaCl()
{
}

void cH2ONaCl:: init_PhaseRegionName()
{
    m_phaseRegion_name[SinglePhase_L]="Single phase(Liquid)";
    m_phaseRegion_name[SinglePhase_V]="Pure vapour phase";
    m_phaseRegion_name[SinglePhase_H]="Pure halite phase";
    m_phaseRegion_name[TwoPhase_V_L]="Vaper + Liquid";
    m_phaseRegion_name[TwoPhase_V_H]="Vaper + Halite";
    m_phaseRegion_name[TwoPhase_L_H]="Liquid + Halite";
    m_phaseRegion_name[ThreePhase_V_L_H]="Vaper + Liquid + Halite";
}
void cH2ONaCl:: init_prop()
{
    m_prop.Region=SinglePhase_L;
    m_prop.T=0;
    m_prop.H=0;
    m_prop.Rho=0;
    m_prop.Rho_l=0;
    m_prop.Rho_v=0;
    m_prop.Rho_h=0;
    m_prop.S_l=0;
    m_prop.S_v=0;
    m_prop.S_h=0;
    m_prop.X_l=0;
    m_prop.X_v=0;
    m_prop.Mu_l=0;
    m_prop.Mu_v=0;
}

void cH2ONaCl:: Calculate()
{
    // 1. 
    m_prop.Region=findRegion(m_T, m_P, m_X);
}

PhaseRegion cH2ONaCl:: findRegion(const double T, const double P, const double X)
{
    double Pres=P/1e5; //Pa -> bar
    double Xl_all=0, Xv_all=0, region_ind=-1;
    // CALCULATE CRITICAL P AND X FOR GIVEN T
    // First we need to find the Critical P and Critical X for the given T
    double cn1[7] = {-2.36, 1.28534e-1, -2.3707e-2, 3.20089e-3, -1.38917e-4, 1.02789e-7, -4.8376e-11};
    double cn2[4] = {2.36, -1.31417e-2, 2.98491e-3, -1.30114e-4};
    double ca[11] = {1, 1.5, 2, 2.5, 3, 4, 5, 1, 2, 2.5, 3};
    double Tcrit_h2o = 373.976;
    double Pcrit_h2o_point = 220.54915;
    double T_der[2] = {499.999, 500};
    double P_der[2] = {0, 0};
    for(size_t i=0;i<2;i++)
    {
        double S=0;
        for(size_t j=0;j<4;j++)
        {
            S+=cn2[j]*pow(T_der[i] - Tcrit_h2o, ca[7+j]);
        }
        P_der[i]=Pcrit_h2o_point + S;
        cout<<P_der[i]<<endl;
    }
    double P_crit = 0;
    double X_crit = 0;  // mole fraction
    double P_crit_h20=0;
    // 2
    if(T<=Tcrit_h2o)
    {
        double Rho_l, Rho_v, h_l, h_v;
        fluidProp_crit_T(T, 1e-8, P_crit_h20,Rho_l, Rho_v, h_l, h_v);
    }
    
    
    return SinglePhase_H;
}

void cH2ONaCl:: fluidProp_crit_T(double T, double tol, double& P,double& Rho_l, double& Rho_v, double& h_l, double& h_v)
{
    T+=273.15; //C -> K

    Rho_l   = 0;
    Rho_v   = 0;
    h_l     = 0;
    h_v     = 0;
    P       = 0;
    double dpd_l=0, dpd_v=0;
    double T_crit = 647.126;
    double T_creg = 646.303775;
    double T_2ph=-1;
    if(T<=T_crit)T_2ph=T;
    bool ind2b=false;
    bool ind2a=false;
    if(T_2ph > T_creg)ind2b=true;
    ind2a=!ind2b;
    if(ind2a)
    {
        int ind2a_iter=1;
        double con_gas = 0.46152200; // kJ/kg
        double T_ind2a = T_2ph;
        double P_ind2a = 0;
        double h_l_ind2a = 0;
        double h_v_ind2a = 0;

        double dpd_l_ind2a = 0;
        double dpd_v_ind2a = 0;
        // [Rho_l_ind2a,Rho_v_ind2a] = approx_Rho_lv(T_ind2a(ind2a_iter));
        double Rho_l_ind2a,Rho_v_ind2a;
        approx_Rho_lv(T_ind2a, Rho_l_ind2a,Rho_v_ind2a);
        cout<<"Rho_l_ind2a: "<<Rho_l_ind2a<<" Rho_v_ind2a: "<<Rho_v_ind2a<<endl;
    }else if(ind2b)
    {

    }else
    {
        cout<<"Never happens in cH2ONaCl:: fluidProp_crit_T"<<endl;
    }
    
}

void cH2ONaCl:: approx_Rho_lv(double T, double& Rho_l , double& Rho_v)
{
    double al1[11] = { 0.3155901296,  -0.03592060636,
                    2.396625841,   -36.39240662,
	                413.6745246,	-2911.342409,
	                12844.66533,	-35543.55367,
	                59925.07856,	-56266.61248,
	                22585.58 };
    double av1[11] = {11.08333753,	-44.64350654,
	                121.0778198,	-278.5153762,
	                327.9464497,	1462.514145,
	                -11593.55916,	38922.19673,
	                -73229.67131,	74466.37149,
	                -31962.35787 };

    double al2[10] = {1.0,	    0.643432247,
	                -25.98811457,	192.8271795,
	                -947.6312526,	3190.638964,
	                -6805.842363,	8088.134131,
	                -4034.574025,	0.0 };
    
    double av2[10] = {1.000000272,  4.026415669,
	                -129.9023268,	1867.667009,
	                -13845.24815,	61824.71587,
	                -171560.6251,	290312.0606,
	                -274292.5181,	111202.097 };
    double xl=0, xv=0;
    double ts = T/647.3;
    double tt=0;
    if(T<=623.15)
    {
        tt = ts - 273.15/647.3;
        cout<<"tt: "<<tt<<endl;
        for(size_t i=0;i<11;i++)
        {
            xl = xl* tt + al1[10-i];
            xv = xv* tt + av1[10-i];
        }
        xv = exp(xv);
        cout<<"xl: "<<xl<<" xv: "<<xv<<endl;
        if(T>623.15)
        {
            tt = pow((1.0 - ts), 0.25);
            for (size_t i = 0; i < 10; i++)
            {
                xl = xl* tt + al2[9-i];
                xv = xv* tt + av2[9-i];
            }
        }
        Rho_l = 1.0/(3.17* xl);
        Rho_v = 1.0/(3.17* xv);
    }
    cout<<"T: "<<T<<endl;
}