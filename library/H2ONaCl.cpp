#include "H2ONaCl.h"

cH2ONaCl::cH2ONaCl(double P, double T, double X)
:m_P(P),
m_T(T),
m_X(Xwt2Xmol(X)),
m_Parray(NULL),
m_Tarray(NULL),
m_Xarray(NULL),
m_number(1),
m_Cr(init_Cr()),
m_f(init_f())
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

f_STRUCT cH2ONaCl:: init_f()
{
    f_STRUCT f={
                {4.64e-3, 5e-7, 1.69078e1, -2.69148e2, 7.63204e3, -4.95636e4, 2.33119e5, -5.13556e5, 5.49708e5, -2.84628e5, 0},
                0
            };
    for(int i=0;i<10;i++)
    {
        f.sum_f10+=f.f[i];
    }
    return f;    
}
Cr_STRUCT cH2ONaCl:: init_Cr()
{
    Cr_STRUCT Cr={
            {
                {0.0, -530.62968529023, 2274.4901424408, -2662.7944829770, 787.79333020687, -69.830527374994},
                {0.0, 17863.832875422, -39514.731563338, 0.0, 33803.884280753, -13855.050202703},
                {6883.3257944332, -256374.36613260, 482125.75981415, 21757.245522644, -341830.16969660, 122231.56417448},
                {0.0, 1179743.3655832, -2173481.0110373, 0.0, 1082995.2168620, -254419.98064049},
                {0.0, -3137777.4947767, 5291191.0757704, -70730.418082074, -1380257.7177877, -251099.14369001},
                {0.0, 4656182.6115608, -7275277.3275387, 0.0, 417742.46148294, 1401635.8244614},
                {0.0, -3155523.1392127, 4792966.6384584, 0.0, 409126.64781209, -1362636.9388386},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 696252.20862664, -1083490.0096447, 0.0, -227228.27401688, 383654.86000660}
            },
            {2, 2, 2, 4},
            {0, 2, 0, 0},
            {-0.225, -1.68, 0.055, -93.0},
            {640, 640, 641.6,270},
            {0.319, 0.319, 0.319, 1.55},
            {34, 40, 30, 1050},
            {20000, 20000, 40000, 25}
    };
    return Cr;
}
void cH2ONaCl:: init_PhaseRegionName()
{
    m_phaseRegion_name[SinglePhase_L]="Single phase(Liquid)";
    m_phaseRegion_name[TwoPhase_L_V_X0]="Liquid + Vapor at X=0";
    m_phaseRegion_name[SinglePhase_V]="Pure vapour phase";
    m_phaseRegion_name[TwoPhase_L_H]="Liquid + Halite";
    m_phaseRegion_name[TwoPhase_V_H]="Vaper + Halite";
    m_phaseRegion_name[ThreePhase_V_L_H]="Vapour + Liquid + Halite";
    m_phaseRegion_name[TwoPhase_V_L_L]="Vapour + Liquid on the liquid side";
    m_phaseRegion_name[TwoPhase_V_L_V]="Vapour + Liquid on the vapour side";
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
    double Xl_all,Xv_all;
    PhaseRegion reg=m_prop.Region=findRegion(m_T, m_P, m_X,Xl_all,Xv_all);
    cout<<m_phaseRegion_name[reg]<<" Xl_all: "<<Xl_all<<" Xv_all: "<<Xv_all<<endl;
}

PhaseRegion cH2ONaCl:: findRegion(const double T, const double P, const double X, double& Xl_all, double& Xv_all)
{
    double Pres=P/1e5; //Pa -> bar
    Xl_all=0, Xv_all=0;
    PhaseRegion region_ind=SinglePhase_L;
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
        // cout<<P_der[i]<<endl;
    }
    
    double P_crit = 0;
    double X_crit = 0;  // mole fraction
    double P_crit_h20=0;
    // 2
    if(T<=Tcrit_h2o)
    {
        double Rho_l, Rho_v, h_l, h_v;
        fluidProp_crit_T(T, 1e-8, P_crit_h20,Rho_l, Rho_v, h_l, h_v);
        P_crit_h20 = P_crit_h20 * 10 ;  // from MPa to Bar
        // this method reproduces Driesner Table of X_v of V+H+L surface, edited by FVehling
        // working vor X_v at V+H+L and at V+H to V  transition
        P_crit=Pcrit_h2o_point + cn1[0]*pow((Tcrit_h2o - T),ca[0]) + cn1[1]*pow((Tcrit_h2o - T),ca[1])
                               + cn1[2]*pow((Tcrit_h2o - T),ca[2]) + cn1[3]*pow((Tcrit_h2o - T),ca[3])
                               + cn1[4]*pow((Tcrit_h2o - T),ca[4]) + cn1[5]*pow((Tcrit_h2o - T),ca[5])
                               + cn1[6]*pow((Tcrit_h2o - T),ca[6]); 
        // cout<<"P_crit: "<<P_crit<<endl;
    }else if(T>Tcrit_h2o && T<=500)
    {
        P_crit = Pcrit_h2o_point + cn2[0]*pow((T - Tcrit_h2o),ca[6+1]) + cn2[1]*pow((T - Tcrit_h2o),ca[6+2])
                                 + cn2[2]*pow((T - Tcrit_h2o),ca[6+3]) + cn2[3]*pow((T - Tcrit_h2o),ca[6+4]);

    }else if(T>500)
    {
        double cn3[3] = {581.0101, (P_der[1]-P_der[0])/(T_der[1]-T_der[0]), -4.88336*1e-4};
        P_crit = cn3[0]*pow((T-500),(11+1 -12)) + cn3[1]*pow((T-500),(11+2 -12))
               + cn3[2]*pow((T-500),(11+3 -12));
    }else
    {
        cout<<"Never happens in cH2ONaCl:: findRegion->P_crit"<<endl;
    }
    
    // X_crit
    double d1[7] = {8e-5, 1e-5, -1.37125e-7, 9.46822e-10, -3.50549e-12, 6.57369e-15, -4.89423e-18};
    double d2[4] = {7.77761e-2, 2.7042e-4, -4.244821e-7, 2.580872e-10};
    if(T>=Tcrit_h2o && T<=600)
    {
        X_crit =  d1[0]*pow((T-Tcrit_h2o),1) + d1[1]*pow((T-Tcrit_h2o),2)
                + d1[2]*pow((T-Tcrit_h2o),3) + d1[3]*pow((T-Tcrit_h2o),4) 
                + d1[4]*pow((T-Tcrit_h2o),5) + d1[5]*pow((T-Tcrit_h2o),6) 
                + d1[6]*pow((T-Tcrit_h2o),7);
    }else if(T>600)
    {
        X_crit = d2[0]*pow((T-600),(1-1)) + d2[1]*pow((T-600),(2-1))
               + d2[2]*pow((T-600),(3-1)) + d2[3]*pow((T-600),(4-1));
    }
    // cout<<"T: "<<T<<" Tcrit_h2o: "<<Tcrit_h2o<<" X_crit: "<<X_crit<<endl;
    // ======================================================================
    // FOR THE V+H & V+L REGION
    double a = 2.4726e-2;
    double b_sub = 1.18061e4;
    double b_boil = 0.941812e4;
    double P_trip_salt = 5e-4;
    double T_trip_salt = 800.7;
    double logP_subboil=0;
    if(T<T_trip_salt)
    {
        logP_subboil= log10(P_trip_salt) + b_sub*(1/(T_trip_salt+273.15) - 1/(T+273.15));
    }else if(T>=T_trip_salt)
    {
        logP_subboil = log10(P_trip_salt) + b_boil*(1/(T_trip_salt+273.15) - 1/(T+273.15));
    }else
    {
        cout<<"Never happens in cH2ONaCl:: findRegion->logP_subboil"<<endl;
    }
    double PNacl = pow(10,(logP_subboil)); // halite vapor pressure
    // cout<<"logP_subboil: "<<logP_subboil<<" PNacl: "<<PNacl<<endl;

    // coeffs
    m_f.f[10] = P_trip_salt - m_f.sum_f10;
    double T_star=0;
    if(T<=T_trip_salt)
    {
        T_star=T/T_trip_salt;
    }
    double P_vlh = m_f.f[0]*(pow(T_star,0)) + m_f.f[1]*(pow(T_star,1)) + m_f.f[2]*(pow(T_star,2)) + m_f.f[3]*(pow(T_star,3))
                 + m_f.f[4]*(pow(T_star,4)) + m_f.f[5]*(pow(T_star,5)) + m_f.f[6]*(pow(T_star,6)) + m_f.f[7]*(pow(T_star,7))
                 + m_f.f[8]*(pow(T_star,8)) + m_f.f[9]*(pow(T_star,9)) + m_f.f[10]*(pow(T_star,10)); 
    if(T>T_trip_salt)//V+L+H suface do not exist
    {
        P_vlh=0; 
    }else if(T==T_trip_salt) //as P_vlh(T_star(T_trip_salt)) is not P_triple_salt
    {
        P_vlh=P_trip_salt;
    }
    // cout<<"T_star: "<<T_star<<" P_vlh:"<<P_vlh<<endl;

    // ======================================================================
    // FOR THE V+L  &  V+H REGION
    // Constants for j-constants
    double k0 = -0.235694;
    double k1 = -0.188838;
    double k2 = 0.004;
    double k3 = 0.0552466;
    double k4 = 0.66918;
    double k5 = 396.848;
    double k6 = 45.0;
    double k7 = -3.2719e-7;
    double k8 = 141.699;
    double k9 = -0.292631;
    double k10 = -0.00139991;
    double k11 = 1.95965e-6;
    double k12 = -7.3653e-10;
    double k13 = 0.904411;
    double k14 = 0.000769766;
    double k15 = -1.18658e-6;
    //  Constants for Xl_vl
    double h1 = 1.68486e-3;
    double h2 = 2.19379e-4;
    double h3 = 4.3858e2;
    double h4 = 1.84508e1;
    double h5 = -5.6765e-10;
    double h6 = 6.73704e-6;
    double h7 = 1.44951e-7;
    double h8 = 3.84904e2;
    double h9 = 7.07477e0;
    double h10 = 6.06896e-5;
    double h11 = 7.62859e-3;
    //  Constants for Xv_vl and Xv_vh
    double j0 = k0 + k1*exp(-k2*T);
    double j1 = k4 + (k3-k4)/(1 + exp((T-k5)/k6)) + k7*(pow((T + k8),2));
    double j2 = k9 + k10*T + k11*(pow(T,2)) + k12*(pow(T,3));
    double j3 = k13 + k14*T + k15*(pow(T,2));
    // cout<<"j0: "<<j0<<" j1: "<<j1<<" j2: "<<j2<<" j3: "<<j3<<endl;
    // ======================================================================
    // Calculate Xl_vlh at V+L+H surface for calculating X_l in V+L region
    double e[6] = {0.0989944 + 3.30796e-6*P_vlh - 4.71759e-10*(pow(P_vlh,2)),
                   0.00947257 - 8.66460e-6*P_vlh + 1.69417e-9*(pow(P_vlh,2)),
                   0.610863 - 1.51716e-5*P_vlh + 1.19290e-8*(pow(P_vlh,2)),        
                   -1.64994 + 2.03441e-4*P_vlh - 6.46015e-8*(pow(P_vlh,2)),
                   3.36474 - 1.54023e-4*P_vlh + 8.17048e-8*(pow(P_vlh,2)),
                   1}; 
    for(int i=0;i<5;i++)e[5]-=e[i];
    // for(int i=0;i<6;i++)cout<<e[i]<<endl;
    double T_hm = T_trip_salt + a*(P_vlh - P_trip_salt);  
    T_star = T/T_hm;
    double Xl_vlh = (e[0]*pow(T_star,0)) + (e[1]*pow(T_star,1)) + (e[2]*pow(T_star,2)) + (e[3]*pow(T_star,3))
            + (e[4]*pow(T_star,4)) + (e[5]*pow(T_star,5));
    if(Xl_vlh>1)Xl_vlh=1;   // X of liquid at the V+L+H surface

    // cout<<"Xl_vlh: "<<Xl_vlh<<endl;

    // ======================================================================
    // Calculate Xl_vh metastable for calulating Xv_vh at V - V+H transition P<P_vlh
    double tol_P_LVH = 1e-6; 
    bool ind=false;
    double Xv_vh=0;
    if(Pres < (P_vlh+tol_P_LVH))
    {
        ind=true;
        // here Pres not P_lvh musst be used? 
        double e2[6] = {0.0989944 + 3.30796e-6*Pres - 4.71759e-10*pow(Pres,2),
                        0.00947257 - 8.66460e-6*Pres + 1.69417e-9*pow(Pres,2),
                        0.610863 - 1.51716e-5*Pres + 1.19290e-8*pow(Pres,2),
                        -1.64994 + 2.03441e-4*Pres - 6.46015e-8*pow(Pres,2),
                        3.36474 - 1.54023e-4*Pres + 8.17048e-8*pow(Pres,2),
                        1};
        for(int i=0;i<5;i++)e2[5]-=e2[i]; 
        double T_hm2 = T_trip_salt + a*(Pres - P_trip_salt);  // here Pres not P_lvh musst be used? 
        double T_star2 = T/T_hm2;
        double Xl_vh = (e2[0]*pow(T_star2,0)) + (e2[1]*pow(T_star2,1)) + (e2[2]*pow(T_star2,2))
                    + (e2[3]*pow(T_star2,3)) + (e2[4]*pow(T_star2,4)) + (e2[5]*pow(T_star2,5));
        // Calculate Xv_vh at V - V+H transition P<P_vlh  
        double P_norm = (Pres - PNacl)/(P_crit - PNacl); // P_crit from equation 5a
        double log10K2 = 1 + j0*(pow((1-P_norm),j1)) + j2*(1-P_norm) + j3*(pow((1-P_norm),2)) - (1+j0+j2+j3)*(pow((1-P_norm),3));  
        double log10K1 = log10K2*(log10(PNacl/P_crit) - log10(Xl_vlh)) + log10(Xl_vlh); // here Xl_vlh must be used, not Xl_vh!?
        double log10K = log10K1 - log10(PNacl/Pres);
        double K_vh = pow(10,log10K);
        Xv_vh = Xl_vh/K_vh;   
    }

    // ======================================================================
    // Calculate Xl_vl in V+L Region 
    double g1 = h2 + (h1-h2)/(1 + exp((T-h3)/h4)) + h5*(T*T);
    double g2 = h7 + (h6-h7)/(1 + exp((T-h8)/h9)) + h10*exp(-h11*T);
    // For Temp<=Tcrit_h2o find X_crit so that for P_crit_h20, which is lower then P_crit for same Temp,  Xl_vl(P_crit_h20) = 0 
    // Note that X_crit is then negative
    // Equation for X_crit comes from eqn for g0 and Xl_vl
    if(T<Tcrit_h2o)
    {
        X_crit =(  ( Xl_vlh - g1*(P_crit - P_vlh) - g2*(pow((P_crit-P_vlh),2)) ) *  sqrt(P_crit-P_crit_h20)/sqrt(P_crit-P_vlh)
          + g1*(P_crit - P_crit_h20) + g2*(pow((P_crit-P_crit_h20),2))  )/( -1 + sqrt(P_crit-P_crit_h20)/(sqrt(P_crit-P_vlh)) ) ;
        // when P_crit < P_crit_h20 for ind_T, then X_crit is complex     
        if(P_crit<P_crit_h20)X_crit=0; 
    }
    // X_crit(P_crit < P_crit_h20) = % this should not happen, but it does near critical point, when P_crit < P_crit_h20  
    // cout<<"X_crit: "<<X_crit<<endl;
    double g0 = (Xl_vlh - X_crit - g1*(P_crit - P_vlh) - g2*pow((P_crit-P_vlh),2))/sqrt(P_crit-P_vlh);
    // cout<<"g0: "<<g0<<endl;

    // if (P_crit < Pres), than Xl_vl is complex 
    double Xl_vl = X_crit + g0*sqrt(P_crit - Pres) + g1*(P_crit - Pres) + g2*(pow((P_crit-Pres),2));  // to low for 1000°C 
    // cout<<"Xl_vl: "<<Xl_vl<<endl;
    //--------------------------------------------------------------------------
    //Calculate Xv_vl in V+L Region  T> T_crit_H2O is ok but constnant minmal
    //offset to Driesner paper
    double P_norm = (Pres - PNacl)/(P_crit - PNacl);
    double log10K2 = 1 + j0*(pow((1-P_norm),j1)) + j2*(1-P_norm) + j3*(pow((1-P_norm),2)) - (1+j0+j2+j3)*(pow((1-P_norm),3));
    double log10K1 = log10K2*(log10(PNacl/P_crit) - log10(Xl_vlh)) + log10(Xl_vlh);
    double log10K = log10K1 - log10(PNacl/Pres);
    double K = pow(10,(log10K));
    double Xv_vl = Xl_vl/K;   // to low mole fraction for 1000°C and 1bar
    if(Pres>P_crit)
    {
        Xv_vl=NAN;
        Xl_vl=NAN;
    }
    if((Pres<=PNacl) && (T>=T_trip_salt))
    {
        Xv_vl=NAN;
        Xl_vl=NAN;
    }
    //--------------------------------------------------------------------------
    //Calculate Regions
    double P_crit_s = P_crit;
    if(P_crit_s < Pcrit_h2o_point)
    {
        P_crit_s=Pcrit_h2o_point; //P=22.141e6 T=375
    }
    // cout<<"P_crit_s: "<<P_crit_s<<endl;
    double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
    double T_crit=0;
    fluidProp_crit_P( Pres*1e5 , 1e-10, T_crit, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8);
    // cout<<"T_crit: "<<T_crit<<endl;
    double Xv = Xv_vl;
    if(ind)Xv = Xv_vh;
    if(Pres>=P_crit)Xv = 0;

    double Xl = Xl_vl;
    // Xl(ind) = 0;
    if(Pres>=P_crit)Xl = 0;
    if(Pres<=P_vlh-tol_P_LVH)Xl = 0;

    if(X<Xv && Pres<=(P_crit_s) && T>=(T_crit) )region_ind  = SinglePhase_V;   // V  & Temp>=(T_crit)
    if(Pres<PNacl && T>T_trip_salt)region_ind = SinglePhase_V;   // V: below NaCl vapor pressure (for all NaCl values)
    if( X == 0 && Pres <= Pcrit_h2o_point && T<(T_crit+1e-9) && T>(T_crit-1e-9))region_ind = TwoPhase_L_V_X0;
    if(X>0 && X>=Xv && T<=T_trip_salt && Pres<=(P_vlh-tol_P_LVH) )region_ind   = TwoPhase_V_H;   // V+H & V+H-surface
    if(X>0 && X>=Xv && T<=T_trip_salt && Pres<(P_vlh+tol_P_LVH) && Pres>(P_vlh-tol_P_LVH) )region_ind   = ThreePhase_V_L_H;
    if(X>0 && X<=Xl && X>=Xv && X>=X_crit && Pres>=(P_vlh+tol_P_LVH)  && Pres<=P_crit_s)region_ind  = TwoPhase_V_L_L;  
    if(X>0 && X>=Xv && X<X_crit && Pres>=(P_vlh+tol_P_LVH) && Pres<=P_crit_s)region_ind           = TwoPhase_V_L_V;  
    //------------------------------------------------------------------------------------------------------
    //FOR THE L+H REGION 
    double ee[6] = {0.0989944 + 3.30796e-6*Pres - 4.71759e-10*(Pres*Pres),
                   0.00947257 - 8.66460e-6*Pres + 1.69417e-9*(Pres*Pres),
                   0.610863 - 1.51716e-5*Pres + 1.19290e-8*(Pres*Pres),
                   -1.64994 + 2.03441e-4*Pres - 6.46015e-8*(Pres*Pres),
                   3.36474 - 1.54023e-4*Pres + 8.17048e-8*(Pres*Pres),
                   1};
    for(int i=0;i<5;i++)ee[5]-=ee[i];
    T_hm = T_trip_salt + a*(Pres - P_trip_salt);  // melting temperature of halite pressure dependent
    double X_hal = (ee[0]*pow((T/T_hm),(1-1))) + (ee[1]*pow((T/T_hm),(2-1))) + (ee[2]*pow((T/T_hm),(3-1)))
    + (e[3]*pow((T/T_hm),(4-1))) + (ee[4]*pow((T/T_hm),(5-1))) + (ee[5]*pow((T/T_hm),(6-1)));

    double X_lh  = X_hal; //Store X of liquid for the L+H region
    if(X>=X_hal &&  T<=T_hm && Pres>=(P_vlh+tol_P_LVH))region_ind = TwoPhase_L_H;  // L+H & L+H-surface
    switch (region_ind)
    {
    case SinglePhase_L:
        Xl_all=X;
        break;
    case TwoPhase_L_H:
        Xl_all=X_lh;
    case ThreePhase_V_L_H:
        Xl_all=Xl;
        break;
    case TwoPhase_V_L_L:
        Xl_all=Xl;
        break;
    case TwoPhase_V_L_V:
        Xl_all=Xl;
        break;
    default:
        break;
    }

    switch (region_ind)
    {
    case SinglePhase_V:
        Xv_all=X;
        break;
    case TwoPhase_V_H:
        Xv_all=Xv;
        break;
    case ThreePhase_V_L_H:
        Xv_all=Xv;
        break;
    case TwoPhase_V_L_L:
        Xv_all=Xv;
        break;
    case TwoPhase_V_L_V:
        Xv_all=Xv;
        break;
    default:
        break;
    }
    return region_ind;
}
void cH2ONaCl:: fluidProp_crit_P(double P, double tol, double& T_2ph, 
                double& Rho_l, double& h_l, double& h_v, double& dpd_l, 
                double& dpd_v, double& Rho_v, double& Mu_l, double& Mu_v)
{
    P = P * 1e-6; //[MPA]
    Rho_l = 0;
    Rho_v = 0;
    h_l = 0;
    h_v = 0;
    T_2ph = 0;
    dpd_l = 0;
    dpd_v = 0;

    double cv = 0;
    double dpT = 0;

    double Tcrit_h2o = 373.976;
    double P_crit = 22.054915;
    double P_creg=21.839129;
    if(P > P_crit )T_2ph=Tcrit_h2o;
    if(P > P_creg && P <= P_crit) //ind2b
    {
        int imax = 20;
        double P_ind2b = P;

        double Rho_l_ind2b = 0;
        double Rho_v_ind2b = 0;
        double h_l_ind2b = 0;
        double h_v_ind2b = 0;
        double dpd_l_ind2b = 0;
        double dpd_v_ind2b = 0;
        double cv_ind2b = 0; 
        double dpT_ind2b = 0;
        //--------------------------------------------------------------------------
        //calculation for 22.05485 > P > P_creg in two phase region
        //--------------------------------------------------------------------------
        double T_ind2b=0;
        if(P_ind2b < 22.05485)
        {
            double T_ind2b_1 = 647.126; 
            double dP = 1;

            double P_ind2b_1 = P_ind2b;

            double Rho_l_ind2b_1 = 0;
            double Rho_v_ind2b_1 = 0;
            double h_l_ind2b_1 = 0;
            double h_v_ind2b_1 = 0;
            double dpd_l_ind2b_1 = 0;
            double dpd_v_ind2b_1 = 0;
            double cv_ind2b_1 = 0;
            double dpT_ind2b_1 = 0;
            double ind2b_iter=1;
            double h_l_c, h_v_c ,Rho_l_c, Rho_v_c, dpd_l_c, dpd_v_c, psa;
            psatc (T_ind2b_1,h_l_c, h_v_c ,Rho_l_c, Rho_v_c, dpd_l_c, dpd_v_c, psa);
            dP = psa - P_ind2b_1;
            Rho_l_ind2b_1 = Rho_l_c;
            Rho_v_ind2b_1 = Rho_v_c;
            h_l_ind2b_1 = h_l_c;
            h_v_ind2b_1 = h_v_c;
            dpd_l_ind2b_1 = dpd_l_c;
            dpd_v_ind2b_1 = dpd_v_c;
            for(int i=0;i<imax;i++)
            {
                double temp, dpsdt;
                approx_ps(T_ind2b_1,temp,dpsdt);
                T_ind2b_1 = T_ind2b_1 - dP/dpsdt;  
                psatc (T_ind2b_1,h_l_c, h_v_c ,Rho_l_c, Rho_v_c, dpd_l_c, dpd_v_c, psa);
                dP = psa - P_ind2b_1;

                Rho_l_ind2b_1 = Rho_l_c;
                Rho_v_ind2b_1 = Rho_v_c;
                h_l_ind2b_1 = h_l_c;
                h_v_ind2b_1 = h_v_c;
                dpd_l_ind2b_1 = dpd_l_c;
                dpd_v_ind2b_1 = dpd_v_c;   
                if(abs(dP) <= tol)
                {
                    break;
                }
            }
            Rho_l_ind2b = Rho_l_ind2b_1;
            Rho_v_ind2b = Rho_v_ind2b_1;

            h_l_ind2b = h_l_ind2b_1;
            h_v_ind2b = h_v_ind2b_1;

            dpd_l_ind2b = dpd_l_ind2b_1;
            dpd_v_ind2b = dpd_v_ind2b_1;

            cv_ind2b = cv_ind2b_1;
            dpT_ind2b = dpT_ind2b_1;

            T_ind2b = T_ind2b_1;
        }else
        {
            double P_ind2b_2 = P_ind2b;

            double Rho_l_ind2b_2 = 0;
            double Rho_v_ind2b_2 = 0;
            double h_l_ind2b_2 = 0;
            double h_v_ind2b_2 = 0;
            double dpd_l_ind2b_2 = 0;
            double dpd_v_ind2b_2 = 0;
            double cv_ind2b_2  = 0;
            double dpT_ind2b_2  = 0;

            double T_ind2b_2 = 0;
            double dP=1;
            double T1 = 647.1259;
            double T2 = 647.126;
            double psa=0;
            for(int i=0;i<imax;i++)
            {
                T_ind2b_2 = 0.5 * (T1 + T2);
                double h_l_c, h_v_c ,Rho_l_c, Rho_v_c, dpd_l_c ,dpd_v_c;
                psatc (T_ind2b_2,h_l_c, h_v_c ,Rho_l_c, Rho_v_c, dpd_l_c ,dpd_v_c, psa);


                if(psa > P_ind2b_2)
                {
                    T2 = T_ind2b_2;
                }else
                {
                    T1 = T_ind2b_2;
                }
                Rho_l_ind2b_2 = Rho_l_c;
                Rho_v_ind2b_2 = Rho_v_c;
                h_l_ind2b_2 = h_l_c;
                h_v_ind2b_2 = h_v_c;

                dpd_l_ind2b_2 = dpd_l_c;
                dpd_v_ind2b_2 = dpd_v_c;
                if(abs(psa - P_ind2b_2) <= tol)break;
            }
            Rho_l_ind2b = Rho_l_ind2b_2;
            Rho_v_ind2b = Rho_v_ind2b_2;

            h_l_ind2b = h_l_ind2b_2;
            h_v_ind2b = h_v_ind2b_2;

            dpd_l_ind2b = dpd_l_ind2b_2;
            dpd_v_ind2b = dpd_v_ind2b_2;

            cv_ind2b = cv_ind2b_2;
            dpT_ind2b = dpT_ind2b_2;

            T_ind2b = T_ind2b_2;
        }
        Rho_l = Rho_l_ind2b;
        Rho_v = Rho_v_ind2b;
        h_l = h_l_ind2b;
        h_v = h_v_ind2b;
        dpd_l = dpd_l_ind2b;
        dpd_v = dpd_v_ind2b;
        cv = cv_ind2b;
        dpT = dpT_ind2b;

        T_2ph = T_ind2b;
    }else if(P<=P_creg) //ind2a
    {
        double P_ind2a=P;
        int imax=20;
        double Pl = 2.302585 + log(P);
        double T_ind2a = 372.83 + Pl * (27.7589 + Pl * (2.3819 + Pl * (0.24834 + Pl * 0.0193855)));
        // cout<<"T_ind2a: "<<T_ind2a<<endl;
        if(T_ind2a<273.15)T_ind2a=273.15;
        if(T_ind2a>647.126)T_ind2a=647.126;
        double psa=0, dpsdt=0;
        approx_ps(T_ind2a, psa, dpsdt);
        // cout<<"psa: "<<psa<<" dpsdt: "<<dpsdt<<endl;exit(0);
        bool ind2a_iter=true;
        if(abs( 1 - psa/P_ind2a) < tol * 10)ind2a_iter=false;

        for(int i=0;i<imax;i++)
        {
            T_ind2a = T_ind2a - (psa - P_ind2a)/dpsdt;  //etwas anders als in water ph
     
            if(T_ind2a < 273.15)
            {
                T_ind2a = 273.15;
            }else if(T_ind2a > 647.126)
            {
                T_ind2a = 647.126;
            }
            approx_ps(T_ind2a,psa, dpsdt);
            if(abs( 1 - psa/P_ind2a) < tol * 10)ind2a_iter=false;
            if(ind2a_iter==false)break;
        }
        // cout<<"T_ind2a: "<<T_ind2a<<"psa: "<<psa<<" dpsdt: "<<dpsdt<<endl;exit(0);
        double al1[11] = {0.3155901296, -0.03592060636, 2.396625841, -36.39240662, 413.6745246, 
                        -2911.342409, 12844.66533, -35543.55367, 59925.07856, -56266.61248, 22585.58};
        double av1[11] = {11.08333753, -44.64350654, 121.0778198, -278.5153762, 327.9464497, 
                        1462.514145, -11593.55916, 38922.19673, -73229.67131, 74466.37149, -31962.35787};
        double al2[10] = {1.0, 0.643432247, -25.98811457, 192.8271795, -947.6312526, 
                        3190.638964, -6805.842363, 8088.134131, -4034.574025, 0.0};
        double av2[10] = {1.000000272, 4.026415669, -129.9023268, 1867.667009, -13845.24815,
                        61824.71587, -171560.6251, 290312.0606, -274292.5181, 111202.097};
        double ts = T_ind2a / 647.3;
        double tt=0;
        double xl = 0, xv = 0;
        if((T_ind2a <= 623.15))
        {
            tt = ts - 273.15 / 647.3;
            
            for(int i=0;i<11;i++)
            {
                xl = xl * tt + al1[10-i];
	            xv = xv * tt + av1[10-i];
            }
            xv = exp(xv);
        }else
        {
            tt = pow((1 - ts), 0.25);
            for(int i=0;i<10;i++)
            {
                xl = xl * tt + al2[9-i];
	            xv = xv * tt + av2[9-i];
            }
        }
        // cout<<"xl: "<<xl<<" xv: "<<xv<<endl;exit(0);
        double Rho_l_ind2a = 1 / (3.17 * xl);
        double Rho_v_ind2a = 1 / (3.17 * xv);

        double h_l_ind2a = 0;
        double h_v_ind2a = 0;
        double dpd_l_ind2a = 0;
        double dpd_v_ind2a = 0;

        double cv_ind2a = 0;
        double dpT_ind2a = 0;
        // %S_l_ind2a = 0;

        double delpl = 1;
        double delpv = 1;
        double delg = 1;

        ind2a_iter=true;
        double con_gas = 0.46152200; //kJ/kg
        for(int i=0;i<imax;i++)
        {
            MP_STRUCT MP = bb(T_ind2a);   
            ID_STRUCT ID = ideal(T_ind2a);
            TWOPHASEPROP_STRUCT l_prop, v_prop;
            twoPhaseProp( T_ind2a,Rho_l_ind2a, Rho_v_ind2a, MP, ID, l_prop,v_prop);
            double dpT_a = (v_prop.s - l_prop.s) * Rho_l_ind2a * Rho_v_ind2a/(Rho_l_ind2a - Rho_v_ind2a) ; 
            // cout<<"dpT_a: "<<dpT_a<<endl;
            delpl = abs(1 - l_prop.p / P_ind2a);
            delpv = abs(1 - v_prop.p / P_ind2a);
            double dt = ( l_prop.f - v_prop.f  + P_ind2a * (1/Rho_l_ind2a - 1/Rho_v_ind2a) )/(l_prop.s - v_prop.s);
            T_ind2a = T_ind2a + dt;
            Rho_l_ind2a = Rho_l_ind2a + (P_ind2a - l_prop.p - l_prop.dpt * dt)/l_prop.dpd;
            Rho_v_ind2a = Rho_v_ind2a + (P_ind2a - v_prop.p - v_prop.dpt * dt)/v_prop.dpd;
            // cout<<"Rho_l_ind2a: "<<Rho_l_ind2a<<endl;
            h_l_ind2a = l_prop.h;
            h_v_ind2a = v_prop.h;
            
            dpd_l_ind2a = l_prop.dpd;
            dpd_v_ind2a = v_prop.dpd;   
            dpT_ind2a = dpT_a; 
            delg = abs( (l_prop.g - v_prop.g)/con_gas/T_ind2a );
            if((delpl < tol) && (delpv < tol) && (delg < tol*1e-2))
            {
                ind2a_iter=false;
                break;
            }
        }
        T_2ph = T_ind2a;

        Rho_l = Rho_l_ind2a;
        Rho_v = Rho_v_ind2a;

        h_l = h_l_ind2a;
        h_v = h_v_ind2a;

        dpd_l = dpd_l_ind2a;
        dpd_v = dpd_v_ind2a;

        cv = cv_ind2a;
        dpT = dpT_ind2a;
        // cout<<"T_2ph: "<<T_2ph<<" Rho_l: "<<Rho_l<<" Rho_v: "<<Rho_v<<" h_l: "<<h_l<<endl;
        // cout<<" h_v: "<<h_v<<" dpd_l: "<<dpd_l<<" dpd_v: "<<dpd_v<<" cv: "<<cv<<" dpT: "<<dpT<<endl;
    }
    //format in si units
    T_2ph = T_2ph -273.15;
    Rho_l = Rho_l * 1e3;
    Rho_v = Rho_v * 1e3;
    h_l = h_l * 1e3;
    h_v = h_v * 1e3;
    dpd_l = dpd_l * 1e3;
    dpd_v = dpd_v * 1e3;
    // cout<<"T_2ph: "<<T_2ph<<" Rho_l: "<<Rho_l<<" Rho_v: "<<Rho_v<<" h_l: "<<h_l<<endl;
    // cout<<" h_v: "<<h_v<<" dpd_l: "<<dpd_l<<" dpd_v: "<<dpd_v<<" cv: "<<cv<<" dpT: "<<dpT<<endl;
}
void cH2ONaCl:: psatc(double T, double& h_l, double& h_v ,double& Rho_l, double& Rho_v,double& dpd_l,double& dpd_v, double& psa)
{
    MP_STRUCT MP = bb(T);   
    ID_STRUCT ID = ideal(T);
    double tt = 1 - (T/647.126);
    double dd = 0.657128 * pow(tt, 0.325);
    double dc = 0.32189;
    Rho_l = dc + dd;
    Rho_v = dc - dd;

    BS_STRUCT BS = base(T, Rho_l, MP);
    RS_STRUCT RS = resid(T, Rho_l);
    TWOPHASEPROP_STRUCT l_prop = props(T,Rho_l,BS,RS,ID);

    BS = base(T, Rho_v, MP);
    RS = resid(T, Rho_v);
    TWOPHASEPROP_STRUCT v_prop = props(T, Rho_v, BS, RS, ID);
    psa = 0.6 * v_prop.p + 0.4 * l_prop.p; 

    h_l = l_prop.h;
    h_v = v_prop.h;

    dpd_l = l_prop.dpd;
    dpd_v = v_prop.dpd;
}
void cH2ONaCl:: approx_ps(double T, double& psa, double& dpsdt)
{
    double a[8] = { -7.8889166, 2.5514255, -6.716169, 33.239495, -105.38479, 174.35319, -148.39348, 48.631602};
    psa = 0;
    dpsdt = 0;
    if(T<=314)
    {
        double pa = 8858.843/T;
        double pb = 607.56335 * pow(T,-0.6);
        psa = 0.1 * exp(6.3573118 - pa + pb);
        dpsdt = psa * (pa - 0.6 * pb)/T;
    }else if(T>314)
    {
        double v = T/647.25;
        double w = abs(1 - v);
        double w2 = sqrt(w);
        double b = 0;
        double c = 0;
        for(int j=0;j<8;j++)
        {
            b = b * w2 + a[7-j];
            c = c * w2 + a[7-j] * 0.5 * (7-j + 2);  //mal nachfragen
        }
        b = b * w/v;
        c = -(c + b)/T;
        psa = 22.093 * exp(b);
        dpsdt = psa * c;
    }
    
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
    //cout<<"ind2a: "<<ind2a<<" ind2b: "<<ind2b<<endl;

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
        double Rho_l_ind2a=0,Rho_v_ind2a=0;
        // cout<<"T_ind2a: "<<T_ind2a<<endl;
        approx_Rho_lv(T_ind2a, Rho_l_ind2a,Rho_v_ind2a);
        // cout<<"Rho_l_ind2a: "<<Rho_l_ind2a<<" Rho_v_ind2a: "<<Rho_v_ind2a<<endl;
        // do while, break when i>20
        // need to optimize
        TWOPHASEPROP_STRUCT l_prop, v_prop;
        for (size_t i = 0; i < 20; i++)
        {
            MP_STRUCT MP = bb(T_ind2a); 
            ID_STRUCT ID = ideal(T_ind2a);
            twoPhaseProp( T_ind2a ,Rho_l_ind2a, Rho_v_ind2a ,MP,ID, l_prop, v_prop);
            
            double delp = abs(1.0 - v_prop.p/l_prop.p);
            double delg = abs((l_prop.g - v_prop.g)/con_gas/T_ind2a);
            double psa = (l_prop.f - v_prop.f)/(1.0/(Rho_v_ind2a) - 1.0/(Rho_l_ind2a));
            //psa could produce higher values than  
            Rho_l_ind2a -= (l_prop.p - psa)/l_prop.dpd;
            Rho_v_ind2a -= (v_prop.p - psa)/v_prop.dpd;
            h_l_ind2a = l_prop.h;
            h_v_ind2a = v_prop.h;
    
            dpd_l_ind2a = l_prop.dpd;
            dpd_v_ind2a = v_prop.dpd;   
    
            P_ind2a = 0.5 * (v_prop.p + l_prop.p);
    
            if((delp < tol) & (delg < tol*1e-2))break;
        }
        Rho_l = Rho_l_ind2a;
        Rho_v = Rho_v_ind2a;
        h_l = h_l_ind2a;
        h_v = h_v_ind2a;

        dpd_l = dpd_l_ind2a; 
        dpd_v = dpd_v_ind2a;

        P = P_ind2a;
        // cout<<"ind2a: Rho_l:"<<Rho_l<<" Rho_v:"<<Rho_v<<" h_l:"
        //     <<h_l<<" h_v:"<<h_v<<" dpd_l:"
        //     <<dpd_l<<" dpd_v:"<<dpd_v<<" P:"<<P<<endl;
    }else if(ind2b)
    {
        double T_ind2b = T_2ph;
  
        MP_STRUCT MP = bb(T_ind2b);   
        ID_STRUCT ID = ideal(T_ind2b);

        double tt = 1.0 - (T_ind2b/647.126);
        double dd = 0.657128 * pow(tt,0.325);
        double dc = 0.32189;
        double Rho_l_ind2b = dc + dd;
        double Rho_v_ind2b = dc - dd;
        TWOPHASEPROP_STRUCT l_prop, v_prop;
        twoPhaseProp( T_ind2b ,Rho_l_ind2b, Rho_v_ind2b ,MP,ID,l_prop,v_prop);

        Rho_l = Rho_l_ind2b;
        Rho_v = Rho_v_ind2b;

        h_l = l_prop.h;
        h_v = v_prop.h;

        dpd_l = l_prop.dpd;
        dpd_v = v_prop.dpd;   

        P = 0.6 * v_prop.p + 0.4 * l_prop.p; //sic
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
        // cout<<"tt: "<<tt<<endl;
        for(size_t i=0;i<11;i++)
        {
            xl = xl* tt + al1[10-i];
            xv = xv* tt + av1[10-i];
        }
        xv = exp(xv);
        // cout<<"xl: "<<xl<<" xv: "<<xv<<endl;
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
    // cout<<"T: "<<T<<endl;
}

MP_STRUCT cH2ONaCl:: bb(double T)
{
    //molecular parameters

    double Cm[2][4] = {{0.7478629, -0.3540782, 0.007159876, -0.003528426},
                        {1.1278334, -0.5944001, -5.010996, 0.63684256}};
    double con_tz = 647.073;  // Kelvin
    double tt = con_tz/T;
    double tt2 = tt*tt;
    double tt3 = tt2 * tt;

    MP_STRUCT MP;
    MP.b1 =  Cm[0][0] + Cm[0][1] * log(1/tt) + (Cm[0][2] + Cm[0][3] * tt2)*tt3;
           
    MP.b1t = (Cm[0][1] - (3* Cm[0][2] + 5* Cm[0][3] * tt2) * tt3 )/T;
        
    MP.b1tt = (-Cm[0][1]  + (12 * Cm[0][2] + 30 * Cm[0][3] * tt2) * tt3 )/T/T;
         
    MP.b2 = Cm[1][0] + ( Cm[1][1] + (Cm[1][2] + Cm[1][3] * tt2) * tt) * tt;
 
    MP.b2t = -(Cm[1][1] + ( 2 * Cm[1][2] + 4 * Cm[1][3] * tt2 ) * tt ) * tt/T;
          
    MP.b2tt = (2 * Cm[1][1]  + (6 * Cm[1][2] + 20 * Cm[1][3] * tt2) * tt) * tt/T/T;
    // cout<<"MP-b1: "<<MP.b1<<" "<<MP.b1t<<" "<<MP.b1tt<<endl;
    // cout<<"MP-b2: "<<MP.b2<<" "<<MP.b2t<<" "<<MP.b2tt<<endl;
    return MP;
}

ID_STRUCT cH2ONaCl::ideal(double T)
{
    //Ideal gas function and its first and (second derivatives)

    double con_gas = 0.46152200;      // kJ/kg
    double con_sref = 7.6180802;      // without Einheit
    double con_uref = -4328.455039;   //K
    double Ci[18] ={19.730271018, 20.9662681977, 
                    -0.483429455355, 6.05743189245, 
                    22.56023885, -9.87532442,
                    -4.3135538513, 0.458155781, 
                    -0.047754901883, 0.0041238460633,
                    -0.00027929052852, 0.000014481695261, 
                    -0.00000056473658748, 0.000000016200446,
                    -0.0000000003303822796, 0.00000000000451916067368,
                    -0.0000000000000370734122708, 0.000000000000000137546068238};

    ID_STRUCT ID;
    double tr = T/100;
    double tl = log(tr);

    ID.f = -con_gas *
        (T * ((Ci[0]/tr + Ci[1]) * tl
        + ((Ci[2]/tr + Ci[3])/tr + Ci[4])/tr + Ci[5]
        + (Ci[6] + (Ci[7] + (Ci[8] + (Ci[9] + (Ci[10] + (Ci[11] 
        + (Ci[12] + (Ci[13] + (Ci[14] + (Ci[15] + (Ci[16]
        + Ci[17] * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr)
            * tr) * tr) * tr + 1 - con_sref) + con_uref);
    
    ID.ft = -con_gas * (Ci[0]/tr + Ci[1] * (tl + 1.0) 
            - (2 * Ci[2]/tr + Ci[3])/tr/tr + Ci[5]
            + (2 * Ci[6] + (3.0 * Ci[7] + (4.0 * Ci[8] + (5.0 * Ci[9]
            + (6 * Ci[10] + (7.0 * Ci[11] + (8.0 * Ci[12]
            + (9 * Ci[13] + (10 * Ci[14] + (11 * Ci[15]
            + (12 * Ci[16] + 13 * Ci[17] * tr) * tr) * tr) * tr) * tr) 
            * tr) * tr) * tr) * tr) * tr) * tr) * tr + 1.0 - con_sref);
    
    ID.ftt = -con_gas/T * (Ci[1] - Ci[0]/tr
	        + (6.0 * Ci[2]/tr + 2.0 * Ci[3])/tr/tr
	        + (2.0 * Ci[6] + (6.0 * Ci[7] + (12.0 * Ci[8]  
            + (20.0 * Ci[9] + (30.0 * Ci[10] + (42.0 * Ci[11]
            + (56.0 * Ci[12] + (72.0 * Ci[13] + (90.0 * Ci[14]
            + (110.0 * Ci[15] + (132.0 * Ci[16] + 156.0 * Ci[17]
            * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr) * tr)
            * tr) * tr);
    // cout<<"ID: "<<ID.f<<" "<<ID.ft<<" "<<ID.ftt<<endl;
    return ID;
}

void cH2ONaCl:: twoPhaseProp(double T,double Rho_l, double Rho_v, MP_STRUCT MP,ID_STRUCT ID, TWOPHASEPROP_STRUCT& l_prop, TWOPHASEPROP_STRUCT& v_prop)
{
    // Computes supposed liquid (mLiq->spro.x) and vapour (mPro->spro.x) properties
    //   and offers difference of gibbs-function delg = |(gl - gv) / (R T)|
    //   bb and ideal must have been called previously


    BS_STRUCT BS = base(T, Rho_l, MP);
    RS_STRUCT RS = resid(T, Rho_l);
    l_prop = props(T,Rho_l,BS,RS,ID);

    BS = base(T, Rho_v, MP);
    RS = resid(T, Rho_v);
    v_prop = props(T,Rho_v,BS,RS,ID);
}

BS_STRUCT cH2ONaCl:: base(double T, double Rho, MP_STRUCT MP)
{
    double con_gas = 0.46152200;   // kJ/kg
    double con_pz  = 0.101325;     // MPa

    double x = 1.0 - 0.25 * MP.b1 * Rho;    //x = saturation
    BS_STRUCT BS;
    BS.x = x;

    double RG_z1 = -log(x) + (91.0 + (-260.0 + 169.0/x)/x)/6.0;
    double RG_z2 = (1.0 + (-130.0 + 169.0/x)/x/3.0)/x;
    double RG_z3 = (1.0 + (-260.0/3.0 + 169.0/x)/x)/x/x;

    double RG_k = MP.b2 - 3.5 * MP.b1;
    double RG_kt = MP.b2t - 3.5 * MP.b1t;
    double RG_ktt = MP.b2tt - 3.5 * MP.b1tt;

    BS.f  = con_gas * T * (RG_z1 + Rho * RG_k + log(Rho * con_gas * T/con_pz));
    BS.fd = con_gas * T * (RG_z2 * MP.b1/4.0 + RG_k + 1.0/Rho);
    BS.fdd = con_gas * T * (RG_z3 * MP.b1 * MP.b1/16.0 - 1.0/Rho/Rho);
    BS.ft = BS.f/T + con_gas * ((T * RG_z2 * MP.b1t / 4.0 + T * RG_kt) * Rho + 1.0);

    BS.ftd = con_gas * (((MP.b1 + T * MP.b1t) * RG_z2
    + T * RG_z3 * MP.b1 * MP.b1t * Rho/4.0)/4.0
    + (RG_k + T * RG_kt) + 1.0/Rho);

    BS.ftt = con_gas  * ((((MP.b1t + T * MP.b1tt / 2) * RG_z2 
	     + T * RG_z3 * MP.b1t * MP.b1t 
	     * Rho / 8)/2 + (2 * RG_kt + T * RG_ktt))
	   * Rho + 1.0/T);
    // cout<<"BS: "<<BS.x<<" "<<BS.f<<" "<<BS.fd<<" "<<BS.fdd<<" "<<BS.ft<<" "<<BS.ftd<<" "<<BS.ftt<<endl;
    return BS;
}

RS_STRUCT cH2ONaCl:: resid(double T, double Rho)
{
    double con_tz = 647.073;  // Kelvin
    // Cr is m_Cr
    double tt = con_tz/T;
    double ed = exp(-Rho);
    double dd = (1 - ed);
    // cout<<"tt: "<<tt<<" ed: "<<ed<<" dd: "<<dd<<endl; exit(0);
    RS_STRUCT RS={0,0,0,0,0,0};
    for(int j=0;j<9;j++)
    {
        double glb_k = m_Cr.g[8-j][0] + 
                        (m_Cr.g[8-j][1] + 
                            (m_Cr.g[8-j][2] + 
                                (m_Cr.g[8-j][3] + 
                                    (m_Cr.g[8-j][4] + m_Cr.g[8-j][5] * tt * tt)* tt
                                ) * tt
                            ) * tt
                        ) * tt;
        // cout<<"glb_k: "<<glb_k<<endl;
        double glb_kt = (m_Cr.g[8-j][1] + 
                            (2 * m_Cr.g[8-j][2] + 
                                (3 * m_Cr.g[8-j][3] + 
                                    (4.0 * m_Cr.g[8-j][4] + 6.0 * m_Cr.g[8-j][5] * tt * tt) * tt
                                ) * tt
                            ) * tt
                        ) * tt;
        // cout<<"glb_kt: "<<glb_kt<<endl;
        double glb_ktt = (2 * m_Cr.g[8-j][1] + 
                            (6 * m_Cr.g[8-j][2] + 
                                (12 * m_Cr.g[8-j][3] + 
                                    (20 * m_Cr.g[8-j][4] + 42 * m_Cr.g[8-j][5] * tt * tt) * tt
                                ) * tt
                            ) * tt
                        ) * tt; 
        // cout<<"glb_ktt: "<<glb_ktt<<endl;
        
        RS.f   = RS.f * dd + glb_k / ((8-j) + 1);
        RS.fd  = RS.fd * dd + glb_k;
        RS.fdd = RS.fdd * dd + glb_k * ((8-j) * ed/dd - 1);
        RS.ft  = RS.ft * dd + glb_kt / ((8-j) + 1);
        RS.ftd = RS.ftd * dd + glb_kt;
        RS.ftt = RS.ftt * dd + glb_ktt/((8-j) + 1); 
        // cout<<"RS.f: "<<RS.f<<" RS.ft: "<<RS.ft<<" RS.ftt: "<<RS.ftt<<endl;
        // cout<<"RS.ftd: "<<RS.ftd<<" RS.fd: "<<RS.fd<<" RS.fdd: "<<RS.fdd<<endl;
        // exit(0);
    }
    // cout<<"RS.f: "<<RS.f<<" RS.ft: "<<RS.ft<<" RS.ftt: "<<RS.ftt<<endl;
    // cout<<"RS.ftd: "<<RS.ftd<<" RS.fd: "<<RS.fd<<" RS.fdd: "<<RS.fdd<<endl;
    // exit(0);
    RS.f   = RS.f * dd;
    RS.fd  = RS.fd * ed;
    RS.fdd = RS.fdd * ed;
    RS.ft  = RS.ft * (-dd/T);
    RS.ftd = RS.ftd * (-ed/T);
    RS.ftt = RS.ftt * dd/T/T;
    // cout<<"RS.f: "<<RS.f<<" RS.ft: "<<RS.ft<<" RS.ftt: "<<RS.ftt<<endl;
    // cout<<"RS.ftd: "<<RS.ftd<<" RS.fd: "<<RS.fd<<" RS.fdd: "<<RS.fdd<<endl;
    // exit(0);
    // local parts
    double tau=0;
    double del=0;
    double dk=0;
    double dl=0;
    for(int j=0;j<4;j++)
    {
        tau = (T - m_Cr.t[j]) / m_Cr.t[j];
        del = (Rho - m_Cr.d[j]) / m_Cr.d[j]; 
        if(abs(del)<1e-9) del= 1e-9; // avoid division by zero

        dk = pow(del, m_Cr.k[j]);
        dl = pow(del, m_Cr.l[j]);
        // cout<<"tau: "<<tau<<" del: "<<del<<" dk: "<<dk<<" dl: "<<dl<<endl;
        double k = m_Cr.gg[j] * dl * exp(-m_Cr.a[j] * dk - m_Cr.b[j] * tau * tau);
        double kt = -2 * m_Cr.b[j] * tau/m_Cr.t[j];
        double kd = (m_Cr.l[j] - m_Cr.a[j] * m_Cr.k[j] * dk)/m_Cr.d[j]/del;
        double kdd = (m_Cr.k[j] * m_Cr.a[j] * dk * (1.0 - m_Cr.k[j]) - m_Cr.l[j])/m_Cr.d[j]/m_Cr.d[j]/del/del;
        double ktt = (4.0 * tau * tau * m_Cr.b[j] - 2.0) * m_Cr.b[j]/m_Cr.t[j]/m_Cr.t[j];    
        // cout<<"k: "<<k<<"kt: "<<kt<<"kd: "<<kd<<" kdd: "<<kdd<<" ktt: "<<ktt<<endl;
        
        RS.f   = RS.f   +  k;
        RS.ft  = RS.ft + ( k * kt );
        RS.ftt = RS.ftt + (k * ktt);
        
        RS.ftd = RS.ftd + (k * kt * kd );
        RS.fd  =  RS.fd  + (k * kd );
        RS.fdd =  RS.fdd + (k * (kdd + kd  * kd) ); 
    }
    // cout<<"RS.f: "<<RS.f<<" RS.ft: "<<RS.ft<<" RS.ftt: "<<RS.ftt<<endl;
    // cout<<"RS.ftd: "<<RS.ftd<<" RS.fd: "<<RS.fd<<" RS.fdd: "<<RS.fdd<<endl;
    return RS;
}
TWOPHASEPROP_STRUCT cH2ONaCl:: props(double T, double Rho, BS_STRUCT BS, RS_STRUCT RS, ID_STRUCT ID)
{
    TWOPHASEPROP_STRUCT PROP;

    PROP.f = BS.f + RS.f + ID.f;
    PROP.p = Rho * Rho * (BS.fd + RS.fd);
    PROP.s = -(BS.ft + RS.ft + ID.ft);

    PROP.g = PROP.f + PROP.p/Rho;
    PROP.u = PROP.f + T * PROP.s;
    PROP.h = PROP.g + T * PROP.s;


    PROP.dpd = 2 * PROP.p/Rho + Rho * Rho * (BS.fdd + RS.fdd);
    PROP.dpt = Rho * Rho * (BS.ftd + RS.ftd);

    PROP.cv = -T * (BS.ftt + RS.ftt + ID.ftt);

    PROP.x = BS.x;
    // cout<<"PROP.f: "<<PROP.f<<" PROP.p: "<<PROP.p<<" PROP.s: "<<PROP.s<<" PROP.g: "<<PROP.g<<" PROP.h: "<<PROP.h<<" PROP.dpd: "<<PROP.dpd<<endl;
    return PROP;
}