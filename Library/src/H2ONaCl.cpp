#include "H2ONaCl.H"
#include "H2ONaCl_LUT_RefineFuncI.H"

// using namespace H2ONaCl;
namespace H2ONaCl
{
    cH2ONaCl::cH2ONaCl()
    :
    // m_Parray(NULL),
    // m_Tarray(NULL),
    // m_Xarray(NULL),
    // m_number(1),
    m_Cr(init_Cr()),
    m_f(init_f()),
    m_colorPrint(false),
    m_num_threads(1),
    m_lut_PTX_2D(NULL),
    m_lut_PTX_3D(NULL)
    {
        // set_num_threads(omp_get_max_threads() > 8 ? 8 : 1);
        init_PhaseRegionName();
        createTable4_Driesner2007a(m_tab4_Driesner2007a);
    }
    cH2ONaCl::~cH2ONaCl()
    {
        destroyLUT();
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
        m_phaseRegion_name[TwoPhase_V_H]="Vapour + Halite";
        m_phaseRegion_name[ThreePhase_V_L_H]="Vapour + Liquid + Halite";
        m_phaseRegion_name[TwoPhase_V_L_L]="Vapour + Liquid on the liquid side";
        m_phaseRegion_name[TwoPhase_V_L_V]="Vapour + Liquid on the vapour side";
        m_phaseRegion_name[UnknownPhaseRegion]="Unknown phase region";

        // m_phaseRegion_name[SinglePhase_L]="0";
        // m_phaseRegion_name[TwoPhase_L_V_X0]="1.5";
        // m_phaseRegion_name[SinglePhase_V]="2";
        // m_phaseRegion_name[TwoPhase_L_H]="3";
        // m_phaseRegion_name[TwoPhase_V_H]="4";
        // m_phaseRegion_name[ThreePhase_V_L_H]="4.5";
        // m_phaseRegion_name[TwoPhase_V_L_L]="5";
        // m_phaseRegion_name[TwoPhase_V_L_V]="6";
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
    void cH2ONaCl::init_prop(H2ONaCl::PROP_H2ONaCl& prop)
    {
        // H2ONaCl::PROP_H2ONaCl prop;
        prop.Region=SinglePhase_L;
        prop.T=0;
        prop.H=0;
        prop.Rho=0;
        prop.Rho_l=0;
        prop.Rho_v=0;
        prop.Rho_h=0;
        prop.S_l=0;
        prop.S_v=0;
        prop.S_h=0;
        prop.X_l=0;
        prop.X_v=0;
        prop.Mu_l=0;
        prop.Mu_v=0;
    }
    ostream & operator << (ostream & out,  cH2ONaCl & A)
    {
        if(A.m_colorPrint)
        {
            out<<COLOR_BLUE<<"--------------------------------------------------------------------------------\n";
            out<<COLOR_RED<<"Temperature = "<<COLOR_DEFAULT<<A.m_prop.T<<" °C, "
               <<COLOR_RED<<"Pressure = "<<COLOR_DEFAULT<<A.m_prop.P/1e5<<" bar, "
               <<COLOR_RED<<"Salinity = "<<COLOR_DEFAULT<<A.m_prop.X_wt*100<<" wt. % NaCl\n";
            out<<COLOR_BLUE<<"--------------------------------------------------------------------------------\n";
            out<<COLOR_GREEN<<"Phase region = "<<COLOR_DEFAULT<<A.m_phaseRegion_name[A.m_prop.Region]<<"\n";
            out<<COLOR_GREEN<<"Bulk density = "<<COLOR_DEFAULT<<A.m_prop.Rho<<" kg/m3, "
               <<COLOR_GREEN<<"Bulk enthalpy = "<<COLOR_DEFAULT<<A.m_prop.H/1000<<" kJ/kg\n";
            out<<COLOR_GREEN<<"Liquid density = "<<COLOR_DEFAULT<<A.m_prop.Rho_l<<" kg/m3\n"
            <<COLOR_GREEN<<"Vapour density = "<<COLOR_DEFAULT<<A.m_prop.Rho_v<<" kg/m3, "
            <<COLOR_GREEN<<"Halite density = "<<COLOR_DEFAULT<<A.m_prop.Rho_h<<" kg/m3\n"
            <<COLOR_GREEN<<"Liquid enthalpy = "<<COLOR_DEFAULT<<A.m_prop.H_l/1000<<" kJ/kg\n"
            <<COLOR_GREEN<<"Vapour enthalpy = "<<COLOR_DEFAULT<<A.m_prop.H_v/1000<<" kJ/kg, "
            <<COLOR_GREEN<<"Halite enthalpy = "<<COLOR_DEFAULT<<A.m_prop.H_h/1000<<" kJ/kg\n"
            <<COLOR_GREEN<<"Liquid saturation = "<<COLOR_DEFAULT<<A.m_prop.S_l<<"\n"
            <<COLOR_GREEN<<"Vapour saturation = "<<COLOR_DEFAULT<<A.m_prop.S_v<<", "
            <<COLOR_GREEN<<"Halite saturation = "<<COLOR_DEFAULT<<A.m_prop.S_h<<"\n"
            <<COLOR_GREEN<<"Liquid viscosity = "<<COLOR_DEFAULT<<A.m_prop.Mu_l<<", "
            <<COLOR_GREEN<<"Vapour viscosity = "<<COLOR_DEFAULT<<A.m_prop.Mu_v<<"\n"
            <<COLOR_GREEN<<"Liquid salinity = "<<COLOR_DEFAULT<<A.m_prop.X_l<<", "
            <<COLOR_GREEN<<"Vapour salinity = "<<COLOR_DEFAULT<<A.m_prop.X_v<<"\n";
            out<<COLOR_BLUE<<"================================================================================\n"<<COLOR_DEFAULT;
        }else
        {
            out<<"--------------------------------------------------------------------------------\n";
            out<<"Temperature = "<<A.m_prop.T<<" °C, "
            <<"Pressure = "<<A.m_prop.P/1e5<<" bar, "
            <<"Salinity = "<<A.m_prop.X_wt*100<<" wt. % NaCl\n";
            out<<"--------------------------------------------------------------------------------\n";
            out<<"Phase region = "<<A.m_phaseRegion_name[A.m_prop.Region]<<"\n";
            out<<"Bulk density = "<<A.m_prop.Rho<<" kg/m3, "
            <<"Bulk enthalpy = "<<A.m_prop.H/1000<<" kJ/kg\n";
            out<<"Liquid density = "<<A.m_prop.Rho_l<<" kg/m3\n"
            <<"Vapour density = "<<A.m_prop.Rho_v<<" kg/m3, "
            <<"Halite density = "<<A.m_prop.Rho_h<<" kg/m3\n"
            <<"Liquid enthalpy = "<<A.m_prop.H_l/1000<<" kJ/kg\n"
            <<"Vapour enthalpy = "<<A.m_prop.H_v/1000<<" kJ/kg, "
            <<"Halite enthalpy = "<<A.m_prop.H_h/1000<<" kJ/kg\n"
            <<"Liquid saturation = "<<A.m_prop.S_l<<"\n"
            <<"Vapour saturation = "<<A.m_prop.S_v<<", "
            <<"Halite saturation = "<<A.m_prop.S_h<<"\n"
            <<"Liquid viscosity = "<<A.m_prop.Mu_l<<", "
            <<"Vapour viscosity = "<<A.m_prop.Mu_v<<"\n"
            <<"Liquid salinity = "<<A.m_prop.X_l<<", "
            <<"Vapour salinity = "<<A.m_prop.X_v<<"\n";
            out<<"================================================================================\n";
            
        }
        
        return out;
    }
    /**
     * \todo 搞清楚pHX的算法流程并补充完剩下的部分
     */
    H2ONaCl::PROP_H2ONaCl cH2ONaCl:: prop_pHX(double p, double H, double X_wt)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        init_prop(prop);
        prop.P=p; prop.H=H; prop.X_wt=X_wt;
        double T1, T2;
        double tol=1e-4;
        guess_T_PhX(p, H, X_wt, T1, T2);
        // cout<<"T1: "<<T1<<" T2: "<<T2<<endl;
        PROP_H2ONaCl prop1, prop2;
        prop1=prop_pTX(p,T1+Kelvin,X_wt, false); 
        prop2=prop_pTX(p,T2+Kelvin,X_wt, false); 
        double h1 = prop1.H;
        double h2 = prop2.H;
        // printf("H: %f X_wt: %f T1: %f T2: %f h1: %f h2: %f\n",H,X_wt,T1,T2, h1, h2);
        double h_l = prop2.H_l;
        PROP_H2ONaCl prop22;
        while(h2 < H)
        {
            T2 = T2 + 1;
            prop22=prop_pTX(p,T2+Kelvin,X_wt, false);
            h2 = prop22.H;
            h_l = prop22.H_l;
            if(h2<H && T2>1000) break;
            if(X_wt==0) break;
        }
        
        if (isnan(T2))
        {
            cout<<"error, T2->prop_pHX is nan: "<<T2<<endl;
            exit(0);
        }
        if(h1>H)
        {
            T1=0;
            prop1=prop_pTX(p,T1+Kelvin,X_wt, false); 
            h1=prop1.H;
        }
        if((T2 > 1000 || h1 > H) || (h2 < H && T2 == 1000  && p >= 1.7e7) )
        {
            prop.Region=UnknownPhaseRegion;
            prop.Rho=NAN;
            prop.Rho_l=NAN;
            prop.Rho_v=NAN;
            prop.Rho_h=NAN;
            prop.H=NAN;
            prop.H_l=NAN;
            prop.H_v=NAN;
            prop.H_h=NAN;
            prop.S_l=NAN;
            prop.S_v=NAN;
            prop.S_h=NAN;
            prop.Mu=NAN;
            prop.Mu_l=NAN;
            prop.X_l=NAN;
            prop.X_v=NAN;
        }
        else
        {
            
            int iteri=0;
            double T_new_mid=0;
            //            double h_search=H;
            bool ind_iter=true;
            //            bool indmu=true;
            while (ind_iter)
            {
                iteri = iteri +1;
                // find new T with regular falsi method
                double T_new=0;
                if(T_new_mid == 0)
                {
                    //                    double dT = (H - h1) * (T2 - T1) / (h2-h1);
                    T_new = T1 +  (H - h1) * (T2 - T1) / (h2-h1); 
                    T_new_mid = 1;
                }else
                {
                    T_new = (T1 +  T2) / 2; 
                    T_new_mid = 0; 
                }
                // claculate new h
                PROP_H2ONaCl PROP_new=prop_pTX(p,T_new+Kelvin,X_wt, false);
                if (isnan(T_new))
                {
                    printf("T_new is nan, T1: %f, T2: %f, H: %f, X:%f, h1: %f, h2: %f\n", T1, T2, H, X_wt, h1, h2);
                    exit(0);
                }
                // cout<<"new H: "<<PROP_new.H<<" region: "<<m_phaseRegion_name[PROP_new.Region]<<endl; exit(0);
                // claculate new h in  L+V+H region 
                calc_sat_lvh(PROP_new, H ,X_wt, false); 
                switch (PROP_new.Region)
                {
                case ThreePhase_V_L_H:
                    {
                        //could happen that S_l is negative, if h is outside of vlh region 
                        PROP_new.H = (PROP_new.S_l * PROP_new.Rho_l * PROP_new.H_l + 
                                    PROP_new.S_v * PROP_new.Rho_v * PROP_new.H_v +
                                    PROP_new.S_h * PROP_new.Rho_h * PROP_new.H_h) / PROP_new.Rho;
                        if(PROP_new.S_l < 0) //calc S_h and S_v 
                        {
                            PROP_new.S_h = (PROP_new.Rho_v * (PROP_new.X_v - X_wt))/(PROP_new.Rho_h * (X_wt-1) + PROP_new.Rho_v * (PROP_new.X_v - X_wt));
                            PROP_new.S_v = 1 - PROP_new.S_h;  
                            PROP_new.Rho = PROP_new.S_v * PROP_new.Rho_v + PROP_new.S_h * PROP_new.Rho_h ;
                            PROP_new.H   = (PROP_new.S_v * PROP_new.Rho_v * PROP_new.H_v + PROP_new.S_h * PROP_new.Rho_h * PROP_new.H_h )/ PROP_new.Rho;
                        }
                        if(PROP_new.S_v<0) //calc S_h and S_l 
                        {
                            PROP_new.S_h = (PROP_new.Rho_l*(PROP_new.X_l-X_wt))/(PROP_new.Rho_h*(X_wt-1) + PROP_new.Rho_l*(PROP_new.X_l-X_wt));
                            PROP_new.S_l = 1 - PROP_new.S_h;  
                            PROP_new.Rho = PROP_new.S_l * PROP_new.Rho_l + PROP_new.S_h * PROP_new.Rho_h ;
                            PROP_new.H   = ( PROP_new.S_l * PROP_new.Rho_l * PROP_new.H_l + PROP_new.S_h * PROP_new.Rho_h * PROP_new.H_h )/ PROP_new.Rho;
                        }
                        if(PROP_new.S_h<0) //calc S_l and S_v 
                        {
                            PROP_new.S_l = (PROP_new.Rho_v*(PROP_new.X_v - X_wt))/(PROP_new.Rho_v*(PROP_new.X_v-X_wt)+ PROP_new.Rho_l*(X_wt-PROP_new.X_l)); 
                            PROP_new.S_v = 1 - PROP_new.S_l;  
                            PROP_new.Rho = PROP_new.S_l * PROP_new.Rho_l + PROP_new.S_v * PROP_new.Rho_v ;
                            PROP_new.H   = ( PROP_new.S_l * PROP_new.Rho_l * PROP_new.H_l + PROP_new.S_v * PROP_new.Rho_v * PROP_new.H_v )/ PROP_new.Rho;                 
                        }
                    }
                    break;
                case TwoPhase_L_V_X0:
                    {
                        // this function has slightly different resutls than fluidprop_TP_Rho
                        // for single pahse X = 0
                        double T_crit, Rho_l, h_l, h_v, dpd_l0, dpd_v0, Rho_v, Mu_l0, Mu_v0;
                        fluidProp_crit_P(p , 1e-12, T_crit, Rho_l, h_l, h_v, dpd_l0, dpd_v0, Rho_v, Mu_l0, Mu_v0);
                        double S_l = (Rho_v*(h_v - H))/(Rho_v*(h_v-H) + Rho_l*(H-h_l)); 
                        double S_v = 1- S_l;
                        double Rho = S_l * Rho_l + S_v * Rho_v; 

                        T_new = T_crit;
                        PROP_new.H     = H;
                        PROP_new.Rho   = Rho;
                        PROP_new.Rho_l = Rho_l;
                        PROP_new.Rho_v = Rho_v;
                        PROP_new.H_l   = h_l;
                        PROP_new.H_v   = h_v;
                        PROP_new.S_l   = S_l;
                        PROP_new.S_v   = S_v;
                        if(S_l>1)
                        {
                            PROP_new.H    = h_l;
                            PROP_new.S_l  = 1;
                            PROP_new.S_v  = 0;
                            PROP_new.Region  = SinglePhase_L;
                            PROP_new.H_v  = 0;
                            PROP_new.Rho_v= 0;
                        }
                        if(S_l<0)
                        {
                            PROP_new.H     = h_v;
                            PROP_new.S_l   = 0;
                            PROP_new.S_v   = 1;
                            PROP_new.Region   = SinglePhase_V;
                            PROP_new.H_l   = 0;
                            PROP_new.Rho_l = 0;
                        }
        
                    }
                    break;
                default:
                    break;
                }
                if(X_wt==1)
                {
                    double X_hal_liq, T_hm;
                    calc_halit_liqidus(p, T_new,X_hal_liq, T_hm);  // T not important 
                    if(T_new <= T_hm && T_new > (T_hm - 1e-4))
                    {
                        double Nenner = ( H * (PROP_new.Rho_l - PROP_new.Rho_h) - (PROP_new.H_l * PROP_new.Rho_l - PROP_new.H_h * PROP_new.Rho_h ) );
                        double S_l_hm = PROP_new.Rho_h * (PROP_new.H_h - H)/ Nenner;
                        double S_h_hm = 1 - S_l_hm;
                        double Rho_hm = S_l_hm * PROP_new.Rho_l + (1-S_l_hm) * PROP_new.Rho_h;
                        double h_hm = ( S_l_hm * PROP_new.Rho_l * PROP_new.H_l + S_h_hm * PROP_new.Rho_h * PROP_new.H_h )/Rho_hm;
                        PROP_new.S_l  = S_l_hm; 
                        PROP_new.S_h  = S_h_hm; 
                        PROP_new.Rho  = Rho_hm; 
                        PROP_new.H    = h_hm; 
                    }
                }
                if (isnan(T_new))
                {
                    cout<<"error, T_new->prop_pHX is nan: "<<T_new<<endl;
                    exit(0);
                }
                
                //  writing global storage
                prop.Region = PROP_new.Region;
                prop.T = T_new;
                prop.H = PROP_new.H;
                prop.Rho = PROP_new.Rho;
                prop.Rho_l = PROP_new.Rho_l;
                prop.Rho_v = PROP_new.Rho_v;
                prop.Rho_h = PROP_new.Rho_h;
                prop.H_l = PROP_new.H_l;
                prop.H_v = PROP_new.H_v;
                prop.H_h = PROP_new.H_h;
                prop.S_l = PROP_new.S_l;
                prop.S_v = PROP_new.S_v;
                prop.S_h = PROP_new.S_h;
                prop.X_l = PROP_new.X_l;
                prop.X_v = PROP_new.X_v;

                if(fabs(PROP_new.H - H)/H  < tol) ind_iter=false;
                if(prop.H> H)
                {
                    T2 = prop.T; 
                    h2= prop.H;  
                }
                if(prop.H< H)
                {
                    T1=prop.T;
                    h1=prop.H;
                }
                // what's this meaning ?????
                //  T1( abs(PROP_new.h - h(ind_iter))./h(ind_iter)  < tol ) = [];  
                // T2( abs(PROP_new.h - h(ind_iter))./h(ind_iter)  < tol ) = [];  
                // h1( abs(PROP_new.h - h(ind_iter))./h(ind_iter)  < tol ) = [];  
                // h2( abs(PROP_new.h - h(ind_iter))./h(ind_iter)  < tol ) = [];   
                
                // ind_iter( abs(PROP_new.h - h(ind_iter))./h(ind_iter)  < tol ) = []; 
                if(iteri>=500)
                {
                    ind_iter=false;
                    break;
                }
            }
            
        }
        if (isnan(prop.T))
        {
            cout<<"error, prop.T is nan in prop_pHX: "<<prop.T<<endl;
            exit(0);
        }
        // calculate dynamic viscosity
        calcViscosity(prop.Region, p, prop.T, prop.X_l, prop.X_v, prop.Mu_l, prop.Mu_v);

        return prop;
    }

    void cH2ONaCl:: calc_halit_liqidus(double Pres, double Temp, double& X_hal_liq, double& T_hm)
    {
        Pres = Pres/1e5;  //[bar]
        double a = 2.4726e-2;
        double P_trip_salt = 5e-4;
        double T_trip_salt = 800.7;

        double e[6] = {0.0989944 + 3.30796e-6*Pres - 4.71759e-10*pow(Pres,2),
                    0.00947257 - 8.66460e-6*Pres + 1.69417e-9*pow(Pres,2),
                    0.610863 - 1.51716e-5*Pres + 1.19290e-8*pow(Pres,2),
                    -1.64994 + 2.03441e-4*Pres - 6.46015e-8*pow(Pres,2),
                    3.36474 - 1.54023e-4*Pres + 8.17048e-8*pow(Pres,2),
                    0
                    };
        e[5]=1 - e[0] - e[1] - e[2] - e[3] - e[4];
        T_hm = T_trip_salt + a*(Pres - P_trip_salt);  // melting temperature of halite is pressure dependent

        X_hal_liq = (e[0]*pow((Temp/T_hm), (1-1))) + 
                    (e[1]*pow((Temp/T_hm), (2-1))) + 
                    (e[2]*pow((Temp/T_hm), (3-1))) +
                    (e[3]*pow((Temp/T_hm), (4-1))) + 
                    (e[4]*pow((Temp/T_hm), (5-1))) + 
                    (e[5]*pow((Temp/T_hm), (6-1)));
    }

    void cH2ONaCl:: calc_sat_lvh(PROP_H2ONaCl& PROP, double h, double X, bool isDeriv)
    {
//        double dh = 1;
//        double dX = 1e6;
        if(PROP.Region==ThreePhase_V_L_H)
        {
            double a = PROP.Rho_l * PROP.X_l - PROP.Rho_l * X;
            double b = PROP.Rho_v * PROP.X_v - PROP.Rho_v * X;
            double c = PROP.Rho_h - PROP.Rho_h * X;

            double e = PROP.Rho_l * PROP.H_l - PROP.Rho_l * h;
            double f = PROP.Rho_v * PROP.H_v - PROP.Rho_v * h;
            double g = PROP.Rho_h * PROP.H_h - PROP.Rho_h * h;

            PROP.S_h = ( a * (f-e) / (b-a)  - e ) /  ( (g - e - ( f-e)* (c-a) / (b-a) ) );
            PROP.S_v = ( - (g-e) * PROP.S_h - e ) / (f-e);
            PROP.S_l = 1 - PROP.S_v - PROP.S_h;
            PROP.Rho =   PROP.S_l * PROP.Rho_l + PROP.S_v * PROP.Rho_v + PROP.S_h * PROP.Rho_h;
            //Deriviation for h
//            if(isDeriv)
//            {
//                double e_plus = PROP.Rho_l * PROP.H_l - PROP.Rho_l * (h +dh);
//                double f_plus = PROP.Rho_v * PROP.H_v - PROP.Rho_v * (h +dh);
//                double g_plus = PROP.Rho_h * PROP.H_h - PROP.Rho_h * (h +dh);
//                double S_h_plus = ( a * (f_plus-e_plus) / (b-a)  - e_plus ) /  ( (g_plus - e_plus - ( f_plus-e_plus)* (c-a) / (b-a) ) );
//                double S_v_plus = ( - (g_plus-e_plus) * S_h_plus - e_plus ) / (f_plus-e_plus); 
//                double S_l_plus = 1 - S_v_plus - S_h_plus;
//            }
        }
    }

    void cH2ONaCl:: guess_T_PhX(double P, double h, double X, double& T1, double& T2)
    {
        P = P*1e-6;
        h = h*1e-3;
        double tol = 1e-6;
        double P_crit = 22.054915;   //MPa
        double P_creg = 21.839129;

        T1=0;
        T2=0;
        
        //1. calculate reg value
        int reg=2; //2 =  two phase pure water region
        if(P>P_crit) reg=4;
        if(reg == 4 && h > 2086) reg=5;

        if(reg==2)
        {
            double hl = 0;
            double hv = 0;
            double hl1[10] = { 0.8733214860, 0.3139052092, 0.09683295123, 0.02789725252, 0.006097758153,
                            0.0009363482053, 0.00009698796863, 0.000006435851441, 0.0000002467695234, 0.000000004155400906
                            };
            double hv1[10] = {1.191055872, -0.237600508, -0.1495666133, -0.05331935637, -0.01282450167,
                            -0.002106161251, -0.0002309312619, -0.00001609217758, -0.000000642634847, -0.00000001117696881
                            };
            double hl2[10] = {0.9506622997, 1.144019695, -8.981135832, 22.05395845, 10.71497685,
                            -183.1514846, 441.2001644, -517.5797252, 310.7525469, -76.7287569
                            };
            double hv2[10] = {0.8378327686, 2.500557338, -14.09731747, 46.9525944977, -85.16478445,
                            70.61754229, 19.13336011, -94.08040182, 75.79109507, -21.14194941
                            };
            double hlr = 0;
            double hvr = 0;
            double P_ind=P;
            if(P_ind < 7) //ind2a
            {
                double pr_ind2a = log(P_ind/22.055);
                // cout<<"pr_ind2a: "<<pr_ind2a<<endl;
                for (size_t i = 0; i < 10; i++)
                {
                    hlr = hlr* pr_ind2a + hl1[9-i];
                    hvr = hvr* pr_ind2a + hv1[9-i];
                }
                hl = hlr * 2086;
                hv = hvr * 2086;
            }else if (P_ind >= 7 && P_ind <= P_creg)//ind2c
            {
                double pr_ind2c = pow(( 1 - P_ind/22.055 ), 0.25);
                // cout<<"pr_ind2c: "<<pr_ind2c<<endl;
                for (size_t i = 0; i < 10; i++)
                {
                    hlr = hlr * pr_ind2c + hl2[9-i];
                    hvr = hvr * pr_ind2c + hv2[9-i];
                }
                // cout<<"hlr: "<<hlr<<endl;
                // cout<<"hvr: "<<hvr<<endl;
                hl = hlr * 2086;
                hv = hvr * 2086;
            }else //ind2b
            {
                hl = 1975;
                hv = 2235; 
            }
            double reg2 = 2;
            if( h < (hl - 80) ) reg2 = 1 ;
            if( h > (hv + 80) ) reg2= 3 ;

            reg = reg2;
        }
        
        // 2. case calculation based on reg value
        switch (reg)
        {
        case 2:
            {
                double h_l0, h_v0, dpd_l0, dpd_v0, Mu_l0, Mu_v0,Rho_l,Rho_v;
                fluidProp_crit_P(P*1e6 , tol, T1, Rho_l, h_l0, h_v0, dpd_l0, dpd_v0, Rho_v, Mu_l0, Mu_v0);
                if(X > 0.6) T1= T1 - 200*X;
                if(X > 0.8) T1= 2;
                if(X <= 0.1) T1 = T1 - 15;
                T2 = T1 + 30 + X*800; //old
                if(X >= 0.4) T2 = T1 + 30 + X*1200;
                if(X < 1e-4) T2= T2 + 15;
                if(T2> 1000) T2 = 1000;
                // cout<<"T1: "<<T1<<" T2: "<<T2<<endl;
            }
            break;
        case 1:
            {
                T1=0;
                double h_l0, h_v0, dpd_l0, dpd_v0, Mu_l0, Mu_v0,Rho_l,Rho_v;
                fluidProp_crit_P(P*1e6 , tol, T2, Rho_l, h_l0, h_v0, dpd_l0, dpd_v0, Rho_v, Mu_l0, Mu_v0);
                T2 = T2 + 1e-7;
                if(X >= 0.4) T2=  T2 + X *200;
                if(X<0.2 && X >= 0.1) T2=  T2 + X *100;
                if(X>0.2 && X <= 0.4) T2=  T2 + X *100;
                if(X>0.4 && X <= 1)  T2=  T2 + X *500;
                if(T2> 1000) T2= 1000;
            }
            break;
        case 3:
            {
                double h_l0, h_v0, dpd_l0, dpd_v0, Mu_l0, Mu_v0,Rho_l,Rho_v;
                fluidProp_crit_P(P*1e6 , tol, T1, Rho_l, h_l0, h_v0, dpd_l0, dpd_v0, Rho_v, Mu_l0, Mu_v0);
                T1 = T1 - 2e-9; 
                if(X >= 0.1  ) T1= T1 + 28;
                // if(X < 1e-3 ) T1= T1; 
                //T1(reg == 3) = T1(reg == 3) - X(reg == 3)*50;
                T2 = 1000;
            }
            break;
        case 4:
            {
                T1 = 0;
                T2 = 450;
                if(X >= 0.3 && X <= 0.5) T2 =  T2 + X *800;
                if(X >= 0.5 && X <= 1) T2 =  T2 + X *1800;
                if(T2> 1000) T2 = 1000;
            }
            break;
        case 5:
            {
                T1 = 350;
                //T1(reg == 5 & X >= 0.3) =  T1(reg == 5 & X >= 0.3) + X(reg == 5 & X >= 0.3) *800;
                T2 = 1000;
            }
            break;
        default:
            break;
        }
    }

    H2ONaCl::PROP_H2ONaCl cH2ONaCl::prop_pTX(double p, double T_K, double X_wt, bool visc_on)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        init_prop(prop);// initialize it first
        prop.P=p; prop.X_wt=X_wt;
        prop.T=T_K-Kelvin;
        //---------------------------------------------------------
        double T=T_K-Kelvin,Xl_all,Xv_all;
        // 1. 
        prop.Region=findRegion(T, p, Xwt2Xmol(X_wt), Xl_all,Xv_all);
        // printf("prop_pTX(p=%.2f bar, T=%E C, X=%E wt)->findRegion: %s\n",p/1E5, T, X_wt, m_phaseRegion_name[prop.Region].c_str());
        // 2. calculate rho
        // still problematic at high T & low P
        double V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out;
        calcRho(prop.Region, T, p, Xl_all, Xv_all, 
                prop.Rho_l, prop.Rho_v, prop.Rho_h, V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out);
        // 3. calculate enthalpy
        calcEnthalpy(prop.Region, T, p, Xl_all, Xv_all, 
                            prop.H_l, prop.H_v, prop.H_h);
        // printf("prop_pTX->calcEnthalpy(MJ/kg): H_l=%.2f, H_v=%.2f, H_h=%.2f, %s\n",prop.H_l,prop.H_v, prop.H_h, m_phaseRegion_name[prop.Region].c_str());
        // 4. 
        double Xw_l = Xl_all * NaCl::MolarMass / (Xl_all * NaCl::MolarMass + (1-Xl_all) * H2O::MolarMass);
        double Xw_v = Xv_all * NaCl::MolarMass / (Xv_all * NaCl::MolarMass + (1-Xv_all) * H2O::MolarMass);

        // 4. calcViscosity
        if(visc_on) calcViscosity(prop.Region, p, T, Xw_l, Xw_v, prop.Mu_l, prop.Mu_v);


        if(prop.Region==SinglePhase_L)prop.S_l=1;
        double Xw=X_wt;
        //  Calculate saturation of liquid in L+V region
        if(prop.Region==TwoPhase_V_L_L | prop.Region==TwoPhase_V_L_V)
        {
            prop.S_l = (prop.Rho_v *(Xw_v - Xw))/(prop.Rho_v*(Xw_v-Xw) + prop.Rho_l *(Xw-Xw_l)); 
        }
        // Calculate saturation of halite in V+H region
        if(prop.Region==TwoPhase_V_H)
        {
            prop.S_h = (prop.Rho_v*(Xw_v-Xw))/(prop.Rho_h*(Xw-1) + prop.Rho_v*(Xw_v-Xw));
        }
        //  Calculate saturation of halite in L+H region  % does not work for X = 1
        if(prop.Region==TwoPhase_L_H)
        {
            prop.S_h = (prop.Rho_l*(Xw_l-Xw))/(prop.Rho_h*(Xw-1) + prop.Rho_l*(Xw_l-Xw));
        }
        
        if(prop.Region==SinglePhase_V) prop.S_v= 1;
        if(prop.Region==TwoPhase_V_L_L || prop.Region==TwoPhase_V_L_V) prop.S_v= 1 - prop.S_l;
        if(prop.Region==TwoPhase_V_H) prop.S_v= 1 - prop.S_h;
        if(prop.Region==TwoPhase_L_H) prop.S_l= 1 - prop.S_h;
        prop.Rho = prop.S_l*prop.Rho_l + prop.S_v*prop.Rho_v + prop.S_h *prop.Rho_h ;
        prop.H = (prop.S_l*prop.Rho_l*prop.H_l + prop.S_v*prop.Rho_v*prop.H_v + prop.S_h * prop.Rho_h * prop.H_h)/prop.Rho;
        // printf("prop_pTX: Rho=%.2f, Rho_l=%.2f, Rho_v=%.2f, Rho_h=%.2f\n",prop.Rho,prop.Rho_l, prop.Rho_v, prop.Rho_h);
        // printf("prop_pTX: S_l=%.2f, S_v=%.2f, S_h=%.2f\n",prop.S_l, prop.S_v, prop.S_h);
        // v+l+h-region: //TODO: why ????
        if(prop.Region==ThreePhase_V_L_H) prop.S_l= NAN; 
        if(prop.Region==ThreePhase_V_L_H) prop.S_v= NAN; 
        if(prop.Region==ThreePhase_V_L_H) prop.S_h= NAN; 
        if(prop.Region==ThreePhase_V_L_H) prop.Rho= NAN; 
        if(prop.Region==ThreePhase_V_L_H) prop.H= NAN; 
        // v+l-region X = 0;
        if(prop.Region==TwoPhase_L_V_X0) prop.S_l= NAN; 
        if(prop.Region==TwoPhase_L_V_X0) prop.S_v= NAN; 
        if(prop.Region==TwoPhase_L_V_X0) prop.S_h= 0; 
        if(prop.Region==TwoPhase_L_V_X0) prop.Rho = NAN; 
        if(prop.Region==TwoPhase_L_V_X0) prop.H= NAN; 
        
        prop.X_l = Xw_l;
        prop.X_v = Xw_v;

        prop.Mu = prop.S_l*prop.Mu_l + prop.S_v*prop.Mu_v; //need to fix later!!! This is not correct, but to test thermophysical model in OpenFoam, use this at this moment
        // v+l+h-region
        if(prop.Region==ThreePhase_V_L_H) prop.Mu= NAN; 
        // v+l-region X = 0;
        if(prop.Region==TwoPhase_L_V_X0) prop.Mu = NAN; 
        
        return prop;
    }
    
    double cH2ONaCl:: rho_pTX(double p, double T_K, double X_wt)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        init_prop(prop);
        prop.T=T_K-Kelvin;
        prop.P=p;
        prop.X_wt=X_wt;
        // 1. 
        double T=T_K-Kelvin,Xl_all,Xv_all;
        prop.Region=findRegion(T, p, Xwt2Xmol(X_wt), Xl_all,Xv_all);
        // 2. calculate rho
        // still problematic at high T & low P
        double V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out;
        calcRho(prop.Region, T, p, Xl_all, Xv_all, 
                prop.Rho_l, prop.Rho_v, prop.Rho_h, V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out);
        // 4. 
        double Xw_l = Xl_all * NaCl::MolarMass / (Xl_all * NaCl::MolarMass + (1-Xl_all) * H2O::MolarMass);
        double Xw_v = Xv_all * NaCl::MolarMass / (Xv_all * NaCl::MolarMass + (1-Xv_all) * H2O::MolarMass);

        if(prop.Region==SinglePhase_L)prop.S_l=1;
        double Xw=X_wt;
        //  Calculate saturation of liquid in L+V region
        if(prop.Region==TwoPhase_V_L_L | prop.Region==TwoPhase_V_L_V)
        {
            prop.S_l = (prop.Rho_v *(Xw_v - Xw))/(prop.Rho_v*(Xw_v-Xw) + prop.Rho_l *(Xw-Xw_l)); 
        }
        // Calculate saturation of halite in V+H region
        if(prop.Region==TwoPhase_V_H)
        {
            prop.S_h = (prop.Rho_v*(Xw_v-Xw))/(prop.Rho_h*(Xw-1) + prop.Rho_v*(Xw_v-Xw));
        }
        //  Calculate saturation of halite in L+H region  % does not work for X = 1
        if(prop.Region==TwoPhase_L_H)
        {
            prop.S_h = (prop.Rho_l*(Xw_l-Xw))/(prop.Rho_h*(Xw-1) + prop.Rho_l*(Xw_l-Xw));
        }
        
        if(prop.Region==SinglePhase_V) prop.S_v= 1;
        if(prop.Region==TwoPhase_V_L_L || prop.Region==TwoPhase_V_L_V) prop.S_v= 1 - prop.S_l;
        if(prop.Region==TwoPhase_V_H) prop.S_v= 1 - prop.S_h;
        if(prop.Region==TwoPhase_L_H) prop.S_l= 1 - prop.S_h;
        prop.Rho = prop.S_l*prop.Rho_l + prop.S_v*prop.Rho_v + prop.S_h *prop.Rho_h ;
        // v+l+h-region
        if(prop.Region==ThreePhase_V_L_H) prop.Rho= NAN; 
        // v+l-region X = 0;
        if(prop.Region==TwoPhase_L_V_X0) prop.Rho = NAN; 
        
        return prop.Rho;
    }

    double cH2ONaCl:: rho_l_pTX(double p, double T_K, double X_wt)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        init_prop(prop);
        prop.T=T_K-Kelvin;
        prop.P=p;
        prop.X_wt=X_wt;
        // 1. 
        double T=T_K-Kelvin,Xl_all,Xv_all;
        prop.Region=findRegion(T, p, Xwt2Xmol(X_wt), Xl_all,Xv_all);
        // 2. calculate rho
        // still problematic at high T & low P
        double V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out;
        calcRho(prop.Region, T, p, Xl_all, Xv_all, 
                prop.Rho_l, prop.Rho_v, prop.Rho_h, V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out);
        
        return prop.Rho_l;
    }

    double cH2ONaCl:: mu_l_pTX(double p, double T_K, double X_wt)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        init_prop(prop);
        prop.T=T_K-Kelvin;
        prop.P=p;
        prop.X_wt=X_wt;
        // 1. 
        double T=T_K-Kelvin,Xl_all,Xv_all;
        prop.Region=findRegion(T, p, Xwt2Xmol(X_wt), Xl_all,Xv_all);
        // 4. 
        double Xw_l = Xl_all * NaCl::MolarMass / (Xl_all * NaCl::MolarMass + (1-Xl_all) * H2O::MolarMass);
        double Xw_v = Xv_all * NaCl::MolarMass / (Xv_all * NaCl::MolarMass + (1-Xv_all) * H2O::MolarMass);

        // 4. calcViscosity
        calcViscosity(prop.Region, p, T, Xw_l, Xw_v, prop.Mu_l, prop.Mu_v);

        return prop.Mu_l;
    }
    double cH2ONaCl:: mu_pTX(double p, double T_K, double X_wt)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        init_prop(prop);
        prop.T=T_K-Kelvin;
        prop.P=p;
        prop.X_wt=X_wt;
        // 1. 
        double T=T_K-Kelvin,Xl_all,Xv_all;
        prop.Region=findRegion(T, p, Xwt2Xmol(X_wt), Xl_all,Xv_all);
        // 2. calculate rho
        // still problematic at high T & low P
        double V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out;
        calcRho(prop.Region, T, p, Xl_all, Xv_all, 
                prop.Rho_l, prop.Rho_v, prop.Rho_h, V_l_out, V_v_out, T_star_l_out, T_star_v_out, n1_v_out, n2_v_out);
        // 4. 
        double Xw_l = Xl_all * NaCl::MolarMass / (Xl_all * NaCl::MolarMass + (1-Xl_all) * H2O::MolarMass);
        double Xw_v = Xv_all * NaCl::MolarMass / (Xv_all * NaCl::MolarMass + (1-Xv_all) * H2O::MolarMass);

        // 4. calcViscosity
        calcViscosity(prop.Region, p, T, Xw_l, Xw_v, prop.Mu_l, prop.Mu_v);


        if(prop.Region==SinglePhase_L)prop.S_l=1;
        double Xw=X_wt;
        //  Calculate saturation of liquid in L+V region
        if(prop.Region==TwoPhase_V_L_L | prop.Region==TwoPhase_V_L_V)
        {
            prop.S_l = (prop.Rho_v *(Xw_v - Xw))/(prop.Rho_v*(Xw_v-Xw) + prop.Rho_l *(Xw-Xw_l)); 
        }
        // Calculate saturation of halite in V+H region
        if(prop.Region==TwoPhase_V_H)
        {
            prop.S_h = (prop.Rho_v*(Xw_v-Xw))/(prop.Rho_h*(Xw-1) + prop.Rho_v*(Xw_v-Xw));
        }
        //  Calculate saturation of halite in L+H region  % does not work for X = 1
        if(prop.Region==TwoPhase_L_H)
        {
            prop.S_h = (prop.Rho_l*(Xw_l-Xw))/(prop.Rho_h*(Xw-1) + prop.Rho_l*(Xw_l-Xw));
        }
        
        if(prop.Region==SinglePhase_V) prop.S_v= 1;
        if(prop.Region==TwoPhase_V_L_L || prop.Region==TwoPhase_V_L_V) prop.S_v= 1 - prop.S_l;
        if(prop.Region==TwoPhase_V_H) prop.S_v= 1 - prop.S_h;
        if(prop.Region==TwoPhase_L_H) prop.S_l= 1 - prop.S_h;
        prop.Mu = prop.S_l*prop.Mu_l + prop.S_v*prop.Mu_v; //need to fix later!!! This is not correct, but to test thermophysical model in OpenFoam, use this at this moment
        // v+l+h-region
        if(prop.Region==ThreePhase_V_L_H) prop.Mu= NAN; 
        // v+l-region X = 0;
        if(prop.Region==TwoPhase_L_V_X0) prop.Mu = NAN; 
        
        return prop.Mu;
    }
    PhaseRegion cH2ONaCl:: findRegion(const double T, const double P, const double X, double& Xl_all, double& Xv_all)
    {
        double Pres=P/1e5; //Pa -> bar
        static_cast<void>(Xl_all=0), Xv_all=0;
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
        // exit(0);
        double P_crit = 0;
        double X_crit = 0;  // mole fraction
        double P_crit_h20=0;
        // 2
        if(T<=Tcrit_h2o)
        {
            double Rho_l, Rho_v, h_l, h_v;
            fluidProp_crit_T(T, 1e-8, P_crit_h20,Rho_l, Rho_v, h_l, h_v);
            // cout<<"T: "<<T<<" P_crit_h20: "<<P_crit_h20<<endl;exit(0);
            P_crit_h20 = P_crit_h20 * 10 ;  // from MPa to Bar
            // this method reproduces Driesner Table of X_v of V+H+L surface, edited by FVehling
            // working vor X_v at V+H+L and at V+H to V  transition
            P_crit=Pcrit_h2o_point + cn1[0]*pow((Tcrit_h2o - T),ca[0]) + cn1[1]*pow((Tcrit_h2o - T),ca[1])
                                + cn1[2]*pow((Tcrit_h2o - T),ca[2]) + cn1[3]*pow((Tcrit_h2o - T),ca[3])
                                + cn1[4]*pow((Tcrit_h2o - T),ca[4]) + cn1[5]*pow((Tcrit_h2o - T),ca[5])
                                + cn1[6]*pow((Tcrit_h2o - T),ca[6]); 
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
            cout<<"Fatal error in cH2ONaCl:: findRegion->P_crit, T: "<<T<<endl;
        }
        // cout<<"T: "<<T<<" P_crit: "<<P_crit<<" P_crit_h20: "<<P_crit_h20<<endl;exit(0);
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
        // cout<<"T: "<<T<<" Tcrit_h2o: "<<Tcrit_h2o<<" X_crit: "<<X_crit<<endl;exit(0);
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
            cout<<"Fatal error in cH2ONaCl:: findRegion->logP_subboil, T: "<<T<<endl;
        }
        double PNacl = pow(10,(logP_subboil)); // halite vapor pressure
        // cout<<"logP_subboil: "<<logP_subboil<<" PNacl: "<<PNacl<<endl;exit(0);

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
        // cout<<"T_star: "<<T_star<<" P_vlh:"<<P_vlh<<endl;exit(0);

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
        // cout<<"j0: "<<j0<<" j1: "<<j1<<" j2: "<<j2<<" j3: "<<j3<<endl;exit(0);
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

        // cout<<"Xl_vlh: "<<Xl_vlh<<endl;exit(0);

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
        // cout<<" Xv_vh: "<<Xv_vh<<endl;exit(0);
        // ======================================================================
        // Calculate Xl_vl in V+L Region 
        double g1 = h2 + (h1-h2)/(1 + exp((T-h3)/h4)) + h5*(T*T);
        double g2 = h7 + (h6-h7)/(1 + exp((T-h8)/h9)) + h10*exp(-h11*T);
        // cout<<"g1: "<<g1<<" g2: "<<g2<<endl;
        // For Temp<=Tcrit_h2o find X_crit so that for P_crit_h20, which is lower then P_crit for same Temp,  Xl_vl(P_crit_h20) = 0 
        // Note that X_crit is then negative
        // Equation for X_crit comes from eqn for g0 and Xl_vl
        // cout<<"X_crit: "<<X_crit<<endl;
        if(T<Tcrit_h2o)
        {
            // when P_crit < P_crit_h20 for ind_T, then X_crit is complex     
            if(P_crit<P_crit_h20)
            {
                X_crit=0; 
            }else
            {
                X_crit =(
                    ( Xl_vlh - g1*(P_crit - P_vlh) - g2*(pow((P_crit-P_vlh),2)) ) *  
                    sqrt(P_crit-P_crit_h20)/sqrt(P_crit-P_vlh) + 
                    g1*(P_crit - P_crit_h20) + 
                    g2*(pow((P_crit-P_crit_h20),2))  
                    )/( -1 + sqrt(P_crit-P_crit_h20)/(sqrt(P_crit-P_vlh)) ) ;
            }
            
        }
        // X_crit(P_crit < P_crit_h20) = % this should not happen, but it does near critical point, when P_crit < P_crit_h20  
        // cout<<"X_crit: "<<X_crit<<endl;exit(0);
        double g0 = (Xl_vlh - X_crit - g1*(P_crit - P_vlh) - g2*pow((P_crit-P_vlh),2))/sqrt(P_crit-P_vlh);
        // cout<<"g0: "<<g0<<endl;

        // if (P_crit < Pres), than Xl_vl is complex. OpenFOAM will crash if calculate sqrt(negative value), IMPORTANT!!!
        double Xl_vl=0,Xv_vl=0;
        if((Pres>P_crit) || ((Pres<=PNacl) && (T>=T_trip_salt)))
        {
            Xv_vl=NAN;
            Xl_vl=NAN;
        }
        else
        {
            Xl_vl = X_crit + g0*sqrt(P_crit - Pres) + g1*(P_crit - Pres) + g2*(pow((P_crit-Pres),2));  // to low for 1000°C 
            
            //Calculate Xv_vl in V+L Region  T> T_crit_H2O is ok but constnant minmal
            //offset to Driesner paper
            double P_norm = (Pres - PNacl)/(P_crit - PNacl);
            double log10K2 = 1 + j0*(pow((1-P_norm),j1)) + j2*(1-P_norm) + j3*(pow((1-P_norm),2)) - (1+j0+j2+j3)*(pow((1-P_norm),3));
            double log10K1 = log10K2*(log10(PNacl/P_crit) - log10(Xl_vlh)) + log10(Xl_vlh);
            double log10K = log10K1 - log10(PNacl/Pres);
            double K = pow(10,(log10K));
            Xv_vl = Xl_vl/K;   // to low mole fraction for 1000°C and 1bar
        }
        // cout<<"Xl_vl: "<<Xl_vl<<" Xv_vl: "<<Xv_vl<<endl;
        //--------------------------------------------------------------------------
        //Calculate Regions
        double P_crit_s = P_crit;
        if(P_crit_s < Pcrit_h2o_point)
        {
            P_crit_s=Pcrit_h2o_point; //P=22.141e6 T=375
        }
        // cout<<"P_crit_s: "<<P_crit_s<<endl;exit(0);
        double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
        double T_crit=0;
        fluidProp_crit_P( Pres*1e5 , 1e-10, T_crit, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8);
        // cout<<"T_crit: "<<T_crit<<endl;exit(0);
        double Xv = Xv_vl;
        if(ind)Xv = Xv_vh;
        if(Pres>=P_crit)Xv = 0;

        double Xl = Xl_vl;
        // Xl(ind) = 0;
        if(Pres>=P_crit)Xl = 0;
        if(Pres<=P_vlh-tol_P_LVH)Xl = 0;
        // cout<<"Xv: "<<Xv<<" Xl: "<<Xl<<endl;
        // printf("Xv: %f\nXl: %f\nP_crit_s: %f\nT_crit: %f\nP_NaCl_vapor: %f\nNaCl::T_Triple: %f\nX_crit: %f\nP_vlh: %f\n\n",
        //         Xv, Xl, P_crit_s, T_crit, PNacl, T_trip_salt, X_crit, P_vlh);

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
        + (ee[3]*pow((T/T_hm),(4-1))) + (ee[4]*pow((T/T_hm),(5-1))) + (ee[5]*pow((T/T_hm),(6-1)));

        double X_lh  = X_hal; //Store X of liquid for the L+H region
        if(X>=X_hal &&  T<=T_hm && Pres>=(P_vlh+tol_P_LVH))region_ind = TwoPhase_L_H;  // L+H & L+H-surface
        // cout<<"Pres: "<<Pres<<" P_vlh: "<<P_vlh<<" tol_P_LVH: "<<tol_P_LVH<<" ddd: "<<(Pres>=(P_vlh+tol_P_LVH))<<endl;
        // cout<<"X: "<<X<<" X_hal: "<<X_hal<<" X-X_hal: "<<X-X_hal<<endl;
        // cout<<"Pres: "<<Pres<<" T_hm: "<<T_hm<<" T: "<<T<<endl;
        // for(int i=0;i<6;i++)cout<<"e["<<i+1<<"]: "<<ee[i]<<endl;
        // cout<<"region_ind: "<<region_ind<<endl;exit(0);
        switch (region_ind)
        {
        case SinglePhase_L:
            Xl_all=X;
            break;
        case TwoPhase_L_H:
            Xl_all=X_lh;
            break;
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
//                double ind2b_iter=1;
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
//                double dP=1;
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
        // cout<<"ind2a: "<<ind2a<<" ind2b: "<<ind2b<<endl;exit(0);

        if(ind2a)
        {
//            int ind2a_iter=1;
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
            // cout<<"Rho_l_ind2a: "<<Rho_l_ind2a<<" Rho_v_ind2a: "<<Rho_v_ind2a<<endl;exit(0);
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
            cout<<"Fatal error in cH2ONaCl:: fluidProp_crit_T"<<endl;
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
        }else
        {
            tt = pow((1.0 - ts), 0.25);
            for (size_t i = 0; i < 10; i++)
            {
                xl = xl* tt + al2[9-i];
                xv = xv* tt + av2[9-i];
            }
        }
        // cout<<"xl: "<<xl<<" xv: "<<xv<<endl;
        Rho_l = 1.0/(3.17* xl);
        Rho_v = 1.0/(3.17* xv);
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
    void cH2ONaCl:: writeProps2xyz(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<H2ONaCl::PROP_H2ONaCl> props, std::string fname, std::string xTitle, std::string yTitle, std::string zTitle, string delimiter)
    {
        cout<<"Writing results to file ..."<<endl;
        if((x.size()*y.size()*z.size())!=props.size())
        {
            cout<<"ERROR: T.size()*P.size()*X.size() != props.size() in writeProps2VTK\nThe VTK file can not be created correctly."<<endl;
            exit(0);
        }
        ofstream fpout(fname);
        if(!fpout)
        {
            cout<<"ERROR: Can not open file: "<<fname<<endl;
            exit(0);
        }
        fpout<<xTitle<<delimiter
             <<yTitle<<delimiter
             <<zTitle<<delimiter
             <<"Temperature(C)"<<delimiter
             <<"Bulk density(kg/m3)"<<delimiter
             <<"Bulk enthalpy(J/kg)"<<delimiter
             <<"Liquid salinity"<<delimiter
             <<"Vapour salinity"<<delimiter
             <<"Liquid density(kg/m3)"<<delimiter
             <<"Vapour density(kg/m3)"<<delimiter
             <<"Halite density(kg/m3)"<<delimiter
             <<"Liquid enthalpy(kJ/kg)"<<delimiter
             <<"Vapour enthalpy(kJ/kg)"<<delimiter
             <<"Halite enthalpy(kJ/kg)"<<delimiter
             <<"Liquid viscosity(Pa s)"<<delimiter
             <<"Vapour viscosity(Pa s)"<<delimiter
             <<"Phase Region(index)"<<delimiter
             <<"Phase Region(name)"<<delimiter
             <<endl;
        for (size_t i = 0; i < x.size(); i++)
        {
            fpout<<x[i]<<delimiter
                 <<y[i]<<delimiter
                 <<z[i]<<delimiter
                 <<props[i].T<<delimiter
                 <<props[i].Rho<<delimiter
                 <<props[i].H/1000.0<<delimiter
                 <<props[i].X_l<<delimiter
                 <<props[i].X_v<<delimiter
                 <<props[i].Rho_l<<delimiter
                 <<props[i].Rho_v<<delimiter
                 <<props[i].Rho_h<<delimiter
                 <<props[i].H_l/1000.0<<delimiter
                 <<props[i].H_v/1000.0<<delimiter
                 <<props[i].H_h/1000.0<<delimiter
                 <<props[i].Mu_l<<delimiter
                 <<props[i].Mu_v<<delimiter
                 <<props[i].Region<<delimiter
                 <<m_phaseRegion_name[props[i].Region]<<delimiter
                 <<endl;
        }
        fpout.close();
    }
    void cH2ONaCl:: writeProps2VTK(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<H2ONaCl::PROP_H2ONaCl> props, std::string fname, bool isNormalize, std::string xTitle, std::string yTitle, std::string zTitle)
    {
        cout<<"Writing results to file ..."<<endl;
        if((x.size()*y.size()*z.size())!=props.size())
        {
            cout<<"ERROR: T.size()*P.size()*X.size() != props.size() in writeProps2VTK\nThe VTK file can not be created correctly."<<endl;
            exit(0);
        }
        ofstream fpout(fname);
        if(!fpout)
        {
            cout<<"ERROR: Can not open file: "<<fname<<endl;
            exit(0);
        }
        string fname_py=fname+".py";
        //  write vtk head
        fpout<<"# vtk DataFile Version 2.0"<<endl;
        fpout<<"Properties of seawater"<<endl;
        fpout<<"ASCII"<<endl;
        fpout<<"DATASET RECTILINEAR_GRID"<<endl;
        fpout<<"DIMENSIONS "<<x.size()<<" "<<y.size()<<" "<<z.size()<<endl;
        double len_x=1, len_y=1, len_z=1;
        double xMAX=max(x);
        double xMIN=min(x);
        double yMAX=max(y);
        double yMIN=min(y);
        double zMAX=max(z);
        double zMIN=min(z);
        double scale_x=1, scale_y=1, scale_z=1;
        len_x=(xMAX==xMIN ? 1: xMAX-xMIN);
        len_y=(yMAX==yMIN ? 1: yMAX-yMIN);
        len_z=(zMAX==zMIN ? 1: zMAX-zMIN);
        if(!isNormalize) // if not normalize the vtk file, write a python script to better visualize result
        {
            scale_y=len_x/len_y;
            scale_z=len_x/len_z;
            ofstream fout_py(fname_py);
            if(!fout_py)
            {
                cout<<"Warning: cannot generate pvPython script for Paraview. "<<fname_py<<endl;
            }else
            {
                fout_py<<"from paraview.simple import *"<<endl;
                fout_py<<"xHvtk = LegacyVTKReader(FileNames=[\'"<<fname<<"\'])"<<endl;
                fout_py<<"renderView1 = GetActiveViewOrCreate(\'RenderView\')"<<endl;
                fout_py<<"xHvtkDisplay = Show(xHvtk, renderView1)"<<endl;
                fout_py<<"xHvtkDisplay.Representation = \'Surface\'"<<endl;
                fout_py<<"renderView1.AxesGrid.Visibility = 1"<<endl;
                fout_py<<"xHvtkDisplay.Scale = ["<<scale_x<<", "<<scale_y<<", "<<scale_z<<"]"<<endl;
                fout_py<<"renderView1.AxesGrid.DataScale = ["<<scale_x<<", "<<scale_y<<", "<<scale_z<<"]"<<endl;
                fout_py<<"renderView1.AxesGrid.DataBoundsInflateFactor = 0"<<endl;
                fout_py<<"renderView1.AxesGrid.XTitle = \'"<<xTitle<<"\'"<<endl;
                fout_py<<"renderView1.AxesGrid.YTitle = \'"<<yTitle<<"\'"<<endl;
                fout_py<<"renderView1.AxesGrid.ZTitle = \'"<<zTitle<<"\'"<<endl;
                if(x.size()>1)fout_py<<"renderView1.AxesGrid.XTitleFontSize = 16"<<endl;
                if(x.size()>1)fout_py<<"renderView1.AxesGrid.XTitleBold = 1"<<endl;
                if(y.size()>1)fout_py<<"renderView1.AxesGrid.YTitleFontSize = 16"<<endl;
                if(y.size()>1)fout_py<<"renderView1.AxesGrid.YTitleBold = 1"<<endl;
                if(z.size()>1)fout_py<<"renderView1.AxesGrid.ZTitleFontSize = 16"<<endl;
                if(z.size()>1)fout_py<<"renderView1.AxesGrid.ZTitleBold = 1"<<endl;
                // set default data source as PhaseRegion
                fout_py<<"#set default data source as PhaseRegion"<<endl;
                fout_py<<"paraview.simple._DisableFirstRenderCameraReset()"<<endl;
                fout_py<<"legacyVTKReader1 = GetActiveSource()"<<endl;
                fout_py<<"renderView1 = GetActiveViewOrCreate('RenderView')"<<endl;
                fout_py<<"legacyVTKReader1Display = GetDisplayProperties(legacyVTKReader1, view=renderView1)"<<endl;
                fout_py<<"ColorBy(legacyVTKReader1Display, ('POINTS', 'PhaseRegion'))"<<endl;
                fout_py<<"legacyVTKReader1Display.RescaleTransferFunctionToDataRange(True, False)"<<endl;
                fout_py<<"legacyVTKReader1Display.SetScalarBarVisibility(renderView1, True)"<<endl;
                fout_py<<"phaseRegionLUT = GetColorTransferFunction('PhaseRegion')\n"<<endl;
                fout_py<<"renderView1.ResetCamera()"<<endl;
                fout_py.close();
            }
        }
        fpout<<"X_COORDINATES "<<x.size()<<" float"<<endl;
        if(isNormalize)
        {
            for(int i=0;i<x.size();i++)fpout<<(x[i]-xMIN)/len_x<<" ";fpout<<endl;
        }else
        {
            for(int i=0;i<x.size();i++)fpout<<x[i]<<" ";fpout<<endl;
        }
        fpout<<"Y_COORDINATES "<<y.size()<<" float"<<endl;
        if(isNormalize)
        {
            for(int i=0;i<y.size();i++)fpout<<(y[i]-yMIN)/len_y<<" ";fpout<<endl;
        }else
        {
            for(int i=0;i<y.size();i++)fpout<<y[i]<<" ";fpout<<endl;
        }
        fpout<<"Z_COORDINATES "<<z.size()<<" float"<<endl;
        if(isNormalize)
        {
            for(int i=0;i<z.size();i++)fpout<<(z[i]-zMIN)/len_z<<" ";fpout<<endl;
        }else
        {
            for(int i=0;i<z.size();i++)fpout<<z[i]<<" ";fpout<<endl;
        }
        fpout<<"POINT_DATA "<<props.size()<<endl;
        // 1. phase region
        fpout<<"SCALARS PhaseRegion int"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Region<<" ";
        }fpout<<endl;
        // temperature
        fpout<<"SCALARS Temperature double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].T<<" ";
        }fpout<<endl;
        // 2. bulk Rho 
        fpout<<"SCALARS Rho double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Rho<<" ";
        }fpout<<endl;
        // 2. bulk H 
        fpout<<"SCALARS H double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].H<<" ";
        }fpout<<endl;
        // 2. Xl 
        fpout<<"SCALARS Xl double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].X_l<<" "; 
        }fpout<<endl;
        // 2. Xv
        fpout<<"SCALARS Xv double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].X_v<<" ";
        }fpout<<endl;
        // 3. rho_l
        fpout<<"SCALARS Rho_l double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Rho_l<<" ";
        }fpout<<endl;
        // 4. Rho_v
        fpout<<"SCALARS Rho_v double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Rho_v<<" ";
        }fpout<<endl;
        // 5. Rho_h
        fpout<<"SCALARS Rho_h double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Rho_h<<" ";
        }fpout<<endl;
        // 6. H_l
        fpout<<"SCALARS H_l double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].H_l<<" ";
        }fpout<<endl;
        // 7. H_v
        fpout<<"SCALARS H_v double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].H_v<<" ";
        }fpout<<endl;
        // 8. H_h
        fpout<<"SCALARS H_h double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].H_h<<" ";
        }fpout<<endl;
        // liquid viscosity
        fpout<<"SCALARS mu_l double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Mu_l<<" ";
        }fpout<<endl;
        // vapour viscosity
        fpout<<"SCALARS mu_v double"<<endl;
        fpout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<props.size();i++)
        {
            fpout<<props[i].Mu_v<<" ";
        }fpout<<endl;
        fpout.close();
    }
    template <typename T>
    T cH2ONaCl:: max(vector<T> data)
    {
        T result=data[0];
        for(int i=0;i<data.size();i++)
        {
            if(data[i]>result)
            {
                result=data[i];
            }
        }
        return result;
    }
    template <typename T>
    T cH2ONaCl:: min(vector<T> data)
    {
        T result=data[0];
        for(int i=0;i<data.size();i++)
        {
            if(data[i]<result)
            {
                result=data[i];
            }
        }
        return result;
    }
    template <typename T> 
    T cH2ONaCl:: sum_array1d(T* a, int n)
    {
        T sum=0;
        for(int i=0;i<n;i++)sum+=a[i];
        return sum;
    };

    void cH2ONaCl:: calcRho(int reg, double T_in, double P_in, double X_l, double X_v, double& Rho_l, double& Rho_v, double& Rho_h, 
                        double& V_l_out, double& V_v_out, double& T_star_l_out, double& T_star_v_out, 
                        double& n1_v_out, double& n2_v_out)
    {
        P_in = P_in/1e5; //Pa to bar
        Rho_l = 0;
        Rho_v = 0;
        Rho_h = 0;

        V_l_out = 0;
        T_star_l_out = 0; 
        V_v_out = 0;
        T_star_v_out = 0; 

        n1_v_out = 0; 
        n2_v_out = 0; 

        const double mass_h2o = 18.015/1e3; 
        const double mass_salt = 58.443/1e3;
        const double P_crit = 220.5491;  //[bar]
        //Fitting parameters to calculate T*
        double n11  = -54.2958 - 45.7623*exp(-9.44785e-4*P_in);
        double n21  = -2.6142 - 0.000239092*P_in; 
        double n22  = 0.0356828 + 4.37235e-6*P_in + 2.0566e-9*pow(P_in,2);
//        double n300 = 7.60664e6/pow((P_in + 472.051),2);
//        double n301 = -50 - 86.1446*exp(-6.21128e-4*P_in);
//        double n302 = 294.318*exp(-5.66735e-3*P_in);
//        double n310 = (-0.0732761*exp(-2.3772e-3*P_in)) - 5.2948e-5*P_in;
//        double n311 = -47.2747 + 24.3653*exp(-1.25533e-3*P_in);
//        double n312 = -0.278529 + 0.00081381*P_in;
        double n20  = 1 - n21*sqrt(n22);
        double n1_1 = 330.47 + 0.942876*sqrt(P_in) + 0.0817193*P_in - 2.47556e-8*pow(P_in,2) + 3.45052e-10*pow(P_in,3);
        double n10  = n1_1;
        double n2_1 = -0.0370751 + 0.00237723*sqrt(P_in) + 5.42049e-5*P_in + 5.84709e-9*pow(P_in,2) - 5.99373e-13*pow(P_in,3);
        double n23  = n2_1 - n20 - n21*sqrt((1+n22));
        double n12  = - n10 - n11;
        // 
        bool ind_lv=(reg==TwoPhase_L_V_X0);
        bool ind_v=(reg==SinglePhase_V || reg==TwoPhase_V_H || reg==ThreePhase_V_L_H || reg==TwoPhase_V_L_L || reg==TwoPhase_V_L_V);
        bool ind_l=(reg==SinglePhase_L || reg==TwoPhase_L_H || reg==ThreePhase_V_L_H || reg==TwoPhase_V_L_L || reg==TwoPhase_V_L_V);
        bool ind_h=(reg==TwoPhase_L_H || reg==TwoPhase_V_H || reg==ThreePhase_V_L_H);
        // cout<<"ind_lv: "<<ind_lv<<" ind_v: "<<ind_v<<" ind_l: "<<ind_l<<" ind_h: "<<ind_h<<endl;
        if(ind_lv)
        {
            double T_2ph0, h_l0, h_v0, dpd_l0, dpd_v0, Mu_l0, Mu_v0;
            fluidProp_crit_P(P_in*1e5, 1e-12, T_2ph0, Rho_l, h_l0, h_v0, dpd_l0, dpd_v0, Rho_v, Mu_l0, Mu_v0);
        }
        if(ind_v)
        {
            double mass_sol_v = mass_h2o*(1-X_v) + mass_salt*X_v;
            double n1_v = n10 + n11*(1-X_v) + n12*pow((1-X_v),2);
            double n2_v = n20 + n21*sqrt(X_v+n22) + n23*X_v;
            double T_star_v = n1_v + n2_v*T_in; // + D_v;  %only for low pres
            double P_star_v = P_in;
            // double Rho_star_l=water_tp_IAPS84(P_star_v*1e5, T_star_v, 50,dRhodP, h, Mu, 1e-9, true); 
            // SteamState S = freesteam_set_pT(P_star_v*1e5, T_star_v+Kelvin);
            // double Rho_star_v=freesteam_rho(S);
            double Rho_star_v=water_rho_pT(P_star_v*1e5, T_star_v+Kelvin);
            bool ind1 = (Rho_star_v > 321.89 && P_star_v <= P_crit); 
            bool ind2 = (std::isnan(Rho_star_v) && P_star_v <= P_crit);
            while (ind1 || ind2)
            {
                double T_2ph, Rho_l, h_l, h_v, dpd_l, dpd_v, Rho_star_v, Mu_l, Mu_v;
                fluidProp_crit_P(P_star_v*1e5,1e-12,T_2ph, Rho_l, h_l, h_v, dpd_l, dpd_v, Rho_star_v, Mu_l, Mu_v);
                ind1 = (Rho_star_v > 321.89 && P_star_v <= P_crit);
                ind2 = (std::isnan(Rho_star_v) && P_star_v <= P_crit); 
            }
            double Vol_v = 1./Rho_star_v;
            double V_v = Vol_v*mass_h2o;

            Rho_v = (mass_sol_v/V_v);
            
            V_v_out = V_v;
            T_star_v_out = T_star_v;
            n1_v_out = n1_v; 
            n2_v_out = n2_v;
        }
        if(ind_l)
        {
            double mass_sol_l = mass_h2o*(1-X_l) + mass_salt*X_l;

            double n1_l = n10 + n11*(1-X_l) + n12*pow((1-X_l),2);
            double n2_l = n20 + n21*sqrt(X_l+n22) + n23*X_l;
            double T_star_l = n1_l + n2_l*T_in; // + D_l; %only for low pres
            double P_star_l = P_in;
            // SteamState S = freesteam_set_pT(P_star_l*1e5, T_star_l+Kelvin);
            // double Rho_star_l=freesteam_rho(S);
            double Rho_star_l=water_rho_pT(P_star_l*1e5, T_star_l+Kelvin);
            // cout<<"P_star_l: "<<P_star_l<<" T_star_l: "<<T_star_l<<" Rho_star_l: "<<Rho_star_l<<endl;
            double Vol = 1/Rho_star_l;
            double V_l = Vol*mass_h2o;
            bool ind_low=( (Rho_star_l < 321.89 || std::isnan(Rho_star_l)) && P_star_l <= P_crit );
            if(ind_low)
            {
                // [ T_crit, Rho_star_crit_l, ~, ~, ~, ~ ] = fluidprop_crit_P( P_star_l(ind_low)*1e5 , 1e-9 );  % P_star_l = P_l
                double T_crit, Rho_star_crit_l, h_l, h_v, dpd_l, dpd_v, Rho_v, Mu_l, Mu_v;
                fluidProp_crit_P(P_star_l*1e5,1e-9,T_crit, Rho_star_crit_l, h_l, h_v, dpd_l, dpd_v, Rho_v, Mu_l, Mu_v);
                double Vol_l_crit = mass_h2o / Rho_star_crit_l;
                // S = freesteam_set_pT(P_star_l*1e5, T_crit-1+Kelvin);
                // double Rho_star_crit_l_minus=freesteam_rho(S);
                double Rho_star_crit_l_minus=water_rho_pT(P_star_l*1e5, T_crit-1+Kelvin);
                double Vol_l_crit_minus = mass_h2o / Rho_star_crit_l_minus;
                double dVol_ldT = (Vol_l_crit - Vol_l_crit_minus)/1;
                double o1 = dVol_ldT ; //- 3*o2 * T_crit.^2;
                double o0 = Vol_l_crit - o1 * T_crit ;//- o2 * T_crit.^3;
                double V_l_low = o0 + o1 *  T_star_l ;// + o2 *  T_star_l(ind_low).^3;
                V_l = V_l_low;
            }
            bool ind_high=((T_in >=600) && (P_in < 390.147) && (X_l > 0.1));
            if(ind_high)
            {
                double P_390 = 390.147;
                double n11_P  = -54.2958 - 45.7623*exp(-9.44785e-4*P_390);
                double n21_P  = -2.6142 - 0.000239092*P_390; 
                double n22_P  = 0.0356828 + 4.37235e-6*P_390 + 2.0566e-9*pow(P_390,2);
                double n20_P  = 1 - n21_P*sqrt(n22_P);
                double n1_1_P = 330.47 + 0.942876*sqrt(P_390) + 0.0817193*P_390 - 2.47556e-8*pow(P_390,2) + 3.45052e-10*pow(P_390,3);
                double n10_P  = n1_1_P;
                double n2_1_P = -0.0370751 + 0.00237723*sqrt(P_390) + 5.42049e-5*P_390 + 5.84709e-9*pow(P_390,2) - 5.99373e-13*pow(P_390,3);
                double n23_P  = n2_1_P - n20_P - n21_P*sqrt((1+n22_P));
                double n12_P  = - n10_P - n11_P;
                double X_l_ind_l = X_l;
                double n1_l_P = n10_P + n11_P*(1-X_l_ind_l) + n12_P*pow((1-X_l_ind_l),2);
                double n2_l_P = n20_P + n21_P*sqrt(X_l_ind_l+n22_P) + n23_P*X_l_ind_l;
                double T_in_ind_l = T_in;
                double T_star_l_P = n1_l_P + n2_l_P*T_in_ind_l;
                double P_400 = 400;
                double n11_P4  = -54.2958 - 45.7623*exp(-9.44785e-4*P_400);
                double n21_P4  = -2.6142 - 0.000239092*P_400; 
                double n22_P4  = 0.0356828 + 4.37235e-6*P_400 + 2.0566e-9*pow(P_400,2);
                double n20_P4  = 1 - n21_P4*sqrt(n22_P4);
                double n1_1_P4 = 330.47 + 0.942876*sqrt(P_400) + 0.0817193*P_400 - 2.47556e-8*pow(P_400,2) + 3.45052e-10*pow(P_400,3);
                double n10_P4  =  n1_1_P4;
                double n2_1_P4 = -0.0370751 + 0.00237723*sqrt(P_400) + 5.42049e-5*P_400 + 5.84709e-9*pow(P_400,2) - 5.99373e-13*pow(P_400,3);
                double n23_P4  = n2_1_P4 - n20_P4 - n21_P4*sqrt((1+n22_P4));
                double n12_P4  = - n10_P4 - n11_P4;
                X_l_ind_l = X_l;
                double n1_l_P4 = n10_P4 + n11_P4*(1-X_l_ind_l) + n12_P4*pow((1-X_l_ind_l),2);
                double n2_l_P4 = n20_P4 + n21_P4*sqrt(X_l_ind_l+n22_P4) + n23_P4*X_l_ind_l;
                T_in_ind_l = T_in;
                double T_star_l_P4 = n1_l_P4 + n2_l_P4*T_in_ind_l;
                double P_1000 = 1000;
                double n11_P1  = -54.2958 - 45.7623*exp(-9.44785e-4*P_1000);
                double n21_P1  = -2.6142 - 0.000239092*P_1000; 
                double n22_P1  = 0.0356828 + 4.37235e-6*P_1000 + 2.0566e-9*pow(P_1000,2);
                double n20_P1  = 1 - n21_P1*sqrt(n22_P);
                double n1_1_P1 = 330.47 + 0.942876*sqrt(P_1000) + 0.0817193*P_1000 - 2.47556e-8*pow(P_1000,2) + 3.45052e-10*pow(P_1000,3);
                double n10_P1  =  n1_1_P1;
                double n2_1_P1 = -0.0370751 + 0.00237723*sqrt(P_1000) + 5.42049e-5*P_1000 + 5.84709e-9*pow(P_1000,2) - 5.99373e-13*pow(P_1000,3);
                double n23_P1  = n2_1_P1 - n20_P1 - n21_P1*sqrt((1+n22_P1));
                double n12_P1  = - n10_P1 - n11_P1;
                X_l_ind_l = X_l;
                double n1_l_P1 = n10_P1 + n11_P1*(1-X_l_ind_l) + n12_P1*pow((1-X_l_ind_l),2);
                double n2_l_P1 = n20_P1 + n21_P1*sqrt(X_l_ind_l+n22_P1) + n23_P1*X_l_ind_l;
                T_in_ind_l = T_in;
                double T_star_l_P1 = n1_l_P1 + n2_l_P1*T_in_ind_l;
                // S = freesteam_set_pT(P_390*1e5, T_star_l_P+Kelvin);
                // double Rho_l_390=freesteam_rho(S);
                double Rho_l_390=water_rho_pT(P_390*1e5, T_star_l_P+Kelvin);
                double Vol_390 = mass_h2o / Rho_l_390; 
                // S = freesteam_set_pT(P_400*1e5, T_star_l_P4+Kelvin);
                // double Rho_l_400=freesteam_rho(S);
                double Rho_l_400=water_rho_pT(P_400*1e5, T_star_l_P4+Kelvin);
                double Vol_400 = mass_h2o / Rho_l_400;   
                // S = freesteam_set_pT(P_1000*1e5, T_star_l_P1+Kelvin); 
                // double Rho_l_1000=freesteam_rho(S);
                double Rho_l_1000=water_rho_pT(P_1000*1e5, T_star_l_P1+Kelvin); 
                double Vol_1000 = mass_h2o / Rho_l_1000; 
                double dVol_dP = (Vol_400 - Vol_390) / (P_400 - P_390);

                double P_610 = P_1000 - P_390;
                double P_1390 = P_1000 + P_390; 
                double o4 = ( - Vol_390 + Vol_1000 - dVol_dP * (P_610) )/ ( - log(P_1390) + log( 2*P_1000 ) - (P_610/P_1390) ) ; 
                double o5 = dVol_dP - o4 / P_1390;
                double o3 = Vol_390 - o4 * log(P_1390) - o5 * P_390;
                double V_l_ind_high = o3 + o4 * log(P_star_l+P_1000) + o5 * P_star_l;
                V_l = V_l_ind_high;
            }
            Rho_l = (mass_sol_l/V_l);
            V_l_out = V_l;
            T_star_l_out = T_star_l; 
            // cout<<"Rho_l: "<<Rho_l<<" V_l_out: "<<V_l_out<<" T_star_l_out: "<<T_star_l_out<<endl;
        }
        if(ind_h)
        {
            double l0  = 2.1704e3;
            double l1  = -2.4599e-1;
            double l2  = -9.5797e-5;
            double l3  = 5.727e-3;
            double l4  = 2.715e-3;
            double l5  = 733.4;
            double l = l3 + l4*exp(T_in/l5);
            double Rho0_h = l0 + l1*(T_in) + l2*pow(T_in,2);
            Rho_h = Rho0_h + l*P_in;
        }
    }

    void cH2ONaCl:: calcEnthalpy(int reg, double T_in, double P_in, double X_l, double X_v,
            double& h_l, double& h_v, double& h_h)
    {
        double P_crit = 220.5491;

        P_in = P_in/1e5;

        h_l = 0;
        h_v = 0;
        h_h = 0;
        double T_star_v_out = 0;
        //Fitting parameters to calculate T*
        double q11  = -32.1724 + 0.0621255*P_in;
        double q21  = -1.69513 - 4.52781e-4*P_in - 6.04279e-8*pow(P_in,2); 
        double q22  = 0.0612567 + 1.88082e-5*P_in;
        double q1_1 = 47.9048 - 9.36994e-3*P_in;
        double q2_1 = 0.241022 + 3.45087e-5*P_in - 4.28356e-9*pow(P_in,2);
        double q12  = -q11 - q1_1;
        double q10  = -q11 - q12;
        double q20  = 1 - q21*sqrt(q22);
        double q23  = q2_1 - q20 - q21*sqrt(1+q22);
        double T_trip_salt = 800.7;
        // P_trip_salt = 5e-4;
        // mass_salt = 58.443/1e3;
        //--------------------------------------------------------------------------
        //FIND ENTHALPY OF VAPOUR
        //find coeff
        bool ind_lv = ( reg==TwoPhase_L_V_X0 );
        bool ind_v  = ( reg == SinglePhase_V || reg == TwoPhase_V_H || reg ==ThreePhase_V_L_H || reg == TwoPhase_V_L_L || reg == TwoPhase_V_L_V ); 
        bool ind_l  = ( reg == SinglePhase_L ||  reg == TwoPhase_L_H || reg ==ThreePhase_V_L_H || reg == TwoPhase_V_L_L || reg == TwoPhase_V_L_V ); 
        bool ind_h  = (reg==TwoPhase_L_H || reg==TwoPhase_V_H || reg==ThreePhase_V_L_H);
        if(ind_lv)
        {
            double T_2ph0, Rho_l0, dpd_l0, dpd_v0, Rho_v0, Mu_l0, Mu_v0;
            fluidProp_crit_P(P_in*1e5, 1e-12,T_2ph0, Rho_l0, h_l, h_v, dpd_l0, dpd_v0, Rho_v0, Mu_l0, Mu_v0);
        }
        if(ind_v)
        {
            double q1_v = q10 + q11*(1-X_v) + q12*pow((1-X_v),2);
            double q2_v = q20 + q21*sqrt(X_v+q22) + q23*X_v;

            double T_star_v = q1_v + q2_v*T_in; 
            double P_star_v = P_in;
            // SteamState S = freesteam_set_pT(P_star_v*1e5, T_star_v+Kelvin);
            // h_v=freesteam_h(S);
            h_v=water_h_pT(P_star_v*1e5, T_star_v+Kelvin);
            bool ind1 = (h_v < 2.086e6 && P_star_v < P_crit);// & P_star_v > 40); 
            bool ind2 = (std::isnan(h_v) && P_star_v < P_crit);// & P_star_v > 40); 
            while (ind1 || ind2)
            {
                double T_2ph0, Rho_l0, h_l0, dpd_l0, dpd_v0, Rho_v0, Mu_l0, Mu_v0;
                fluidProp_crit_P(P_star_v*1e5, 1e-12,T_2ph0, Rho_l0, h_l0, h_v, dpd_l0, dpd_v0, Rho_v0, Mu_l0, Mu_v0);
                ind1 = (h_v < 2.086e6 && P_star_v < P_crit && P_star_v > 40);
                ind2 = (std::isnan(h_v) && P_star_v < P_crit && P_star_v > 40); 
            }
            T_star_v_out= T_star_v; 
        }
        // FIND ENTHALPY OF LIQUID
        if(ind_l)
        {
            double q1_l = q10 + q11*(1-X_l) + q12*pow((1-X_l),2);
            double q2_l = q20 + q21*sqrt(X_l+q22) + q23*X_l;
            //from Driesner is equal to above formulation
//            double q1_lb = q1_1 + q11 * (1-X_l) - (q1_1 +q11) * pow((1-X_l),2);
//            double q2_lb = 1 - q21 * sqrt(q22) + q21 * sqrt(X_l+q22) + X_l * (q21 * sqrt(q22) - 1- q21 * sqrt(1+q22) + q2_1);
            double T_star_l = q1_l + q2_l*T_in;
            double P_star_l = P_in;
            // SteamState S = freesteam_set_pT(P_star_l*1e5, T_star_l+Kelvin);
            // h_l=freesteam_h(S);
            h_l=water_h_pT(P_star_l*1e5, T_star_l+Kelvin);
            // printf("P_star_l: %f, T_star_l: %f, h_l: %f\n",P_star_l, T_star_l,h_l);
            //nedded for boiling temps from 180 to Tcri, is not in Driesners Paper
            bool ind_low = ( (h_l > 2.086e6 || std::isnan(h_l))  &&  P_star_l < P_crit  &&  T_in < 375 );
            if(ind_low)
            {
                //  find boiling temperature and spec enthlapy there for given Pressure
                double T_crit, Rho_l0, h_l_crit, dpd_l0, dpd_v0, Rho_v0, Mu_l0, Mu_v0;
                fluidProp_crit_P(P_star_l*1e5, 1e-9,T_crit, Rho_l0, h_l_crit, h_v, dpd_l0, dpd_v0, Rho_v0, Mu_l0, Mu_v0);//P_star_l = P_l
                // find derivative of spec enthlapy at boiling temperature for given Pressure
                // S = freesteam_set_pT(P_star_l*1e5, T_crit-1+Kelvin);
                // double h_l_minus=freesteam_h(S);
                double h_l_minus=water_h_pT(P_star_l*1e5, T_crit-1+Kelvin);
                double dh_ldT = (h_l_crit - h_l_minus)/1;
                double o1 = dh_ldT ; 
                double o0 = h_l_crit - o1 * T_crit ;
                h_l = o0 + o1 *  T_star_l;
            }
            // printf("h_l: %f\n",h_l);
            bool ind_high = ( P_in <= 390.147  &&  T_in > 600);
            if(ind_high)
            {
                double P_390 = 390.147;
                // const for P = 390.147 bar
                q11  = -7.934322551500003;
                q21  = -1.880979162365801; 
                q22  =  0.068594662805400;
                q12  = -36.31482346732;
                q10  =  44.24914601882;
                q20  =  1.492639405203696;
                q23  =  0.705615854382021;
                q1_l = q10 + q11*(1-X_l) + q12*pow((1-X_l),2);
                q2_l = q20 + q21*sqrt(X_l+q22) + q23*X_l;
                double T_star_l_P390 = q1_l + q2_l*T_in;

                double P4 = 400;
                q11  = -7.322200000000002;
                q21  = -1.885910864000000; 
                q22  = 0.068779980000000;
                q12  = -36.834624000000000;
                q10  = 44.156824000000000;
                q20  = 1.494597805305935;
                q23  = 0.709231197191129;
                q1_l = q10 + q11*(1-X_l) + q12*pow((1-X_l),2);
                q2_l = q20 + q21*sqrt(X_l+q22) + q23*X_l;
                double T_star_l_P4 = q1_l + q2_l*T_in;
                
                double P1 = 1000;
                q11  = 29.953100000000000;
                q21  = -2.208338900000000; 
                q22  = 0.080064900000000;
                q12  = -68.487960000000000;
                q10  = 38.534860000000000;
                q20  = 1.624865871647275;
                q23  = 0.941423327837196;
                q1_l = q10 + q11*(1-X_l) + q12*pow((1-X_l),2);
                q2_l = q20 + q21*sqrt(X_l+q22) + q23*X_l;
               double T_star_l_P1 = q1_l + q2_l*T_in;

                // S = freesteam_set_pT(P_390*1e5, T_star_l_P390+Kelvin);
                // double h_l_390=freesteam_h(S);
                double h_l_390=water_h_pT(P_390*1e5, T_star_l_P390+Kelvin);

                // S = freesteam_set_pT(P4*1e5, T_star_l_P4+Kelvin);
                // double h_l_400=freesteam_h(S);
                double h_l_400=water_h_pT(P4*1e5, T_star_l_P4+Kelvin);

                // S = freesteam_set_pT(P1*1e5, P1+Kelvin);
                // double h_l_1000=freesteam_h(S);
                double h_l_1000=water_h_pT(P1*1e5, T_star_l_P1+Kelvin);
                // printf("hl390: %f, hl400: %f, hl1000: %f\n",h_l_390, h_l_400, h_l_1000);
                double dh_l_dP = (h_l_400 - h_l_390) / (P4 - P_390);
                double P_610 = P1 - P_390;
                double P_1390 = P1 + P_390; 
                double o4 = ( - h_l_390 + h_l_1000 - dh_l_dP * (P_610) )/( - log(P_1390) + log( 2*P1 ) - (P_610/P_1390) ) ; 
                double o5 = dh_l_dP - o4 / P_1390;
                double o3 = h_l_390 - o4 * log(P_1390) - o5 * P_390;
                double h_l_ind_high = o3 + o4 * log(P_star_l+P1) + o5 * P_star_l;
                // printf("h_l_ind_high: %f, o3: %f, o4: %f, P_star_l: %f, P1: %f, o5: %f\n",h_l_ind_high, o3, o4, P_star_l, P1, o5);
                h_l = h_l_ind_high;
            }
        }
        if(ind_h)
        {
            h_h = m_NaCl.SpecificEnthalpy(T_in,P_in);
        }
    }

    void cH2ONaCl:: calcViscosity(int reg, double P, double T, double Xw_l, double Xw_v, double& mu_l, double& mu_v)
    {
        double a1 = -35.9858;
        double a2 = 0.80017;
        double b1 = 1e-6;
        double b2 = -0.05239;
        double b3 = 1.32936;
        mu_l = 0;
        mu_v = 0;
        // calculation of mu liquid
        bool ind_l=(reg==SinglePhase_L || reg==TwoPhase_L_V_X0 || reg==TwoPhase_L_H || reg==ThreePhase_V_L_H || reg==TwoPhase_V_L_L || reg==TwoPhase_V_L_V);
        if(ind_l)
        {
            double e1 = a1 * pow(Xw_l,a2);
            double e2 = 1 - b1 * pow(T,b2) - b3 * pow(Xw_l,a2) * pow(T,b2); 
            double T_star_l = e1 + e2 * T;
            if(std::isnan(T_star_l))T_star_l = 0;
            // SteamState S = freesteam_set_pT(P, T_star_l+Kelvin);
            // mu_l=freesteam_mu(S);
            mu_l=water_mu_pT(P, T_star_l+Kelvin);
            if(std::isnan(mu_l))
            {
                double T_2ph0, Rho_l0, h_l0,h_v0, dpd_l0, dpd_v0, Rho_v0, Mu_v0;
                fluidProp_crit_P(P, 1e-10,T_2ph0, Rho_l0, h_l0, h_v0, dpd_l0, dpd_v0, Rho_v0, mu_l, Mu_v0);
            }
        }
        bool ind_v = ( reg==TwoPhase_L_V_X0 | reg==SinglePhase_V | reg==TwoPhase_V_H | reg==ThreePhase_V_L_H | reg==TwoPhase_V_L_L | reg==TwoPhase_V_L_V);
        if(ind_v)
        {
            double e1 = a1 * pow(Xw_v,a2);
            double e2 = 1 - b1 * pow(T,b2) - b3 * pow(Xw_v,a2) * pow(T,b2); 
            double T_star_v = e1 + e2 * T;
            
            bool ind_0 = (T_star_v > 0);
            if(ind_0)
            {
                // SteamState S = freesteam_set_pT(P, T_star_v+Kelvin);
                // mu_v=freesteam_mu(S);
                mu_v=water_mu_pT(P, T_star_v+Kelvin);
            }
            if(std::isnan(mu_v))
            {
                double T_2ph0, Rho_l0, h_l0,h_v0, dpd_l0, dpd_v0, Rho_v0, mu_l0;
                fluidProp_crit_P(P, 1e-10,T_2ph0, Rho_l0, h_l0, h_v0, dpd_l0, dpd_v0, Rho_v0, mu_l0, mu_v);
            }
        }
    }
    double cH2ONaCl::water_rho_pT(double p, double T_K)
    {
        #ifdef USE_PROST
            double d, dp, ds, dh;
            Prop *prop0;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            prop0 = newProp('t', 'p', 1);
            d = 0.0;
            water_tp(T_K,p,d,dp,prop0);
            d=prop0->d;
            // very very important!!!!
            prop0 = freeProp(prop0);
            return d;
        #else 
            SteamState S = freesteam_set_pT(p, T_K);
            return freesteam_rho(S);
        #endif 
    }

    double cH2ONaCl::water_h_pT(double p, double T_K)
    {
        #ifdef USE_PROST
            double d, dp, ds, dh;
            Prop *prop0;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            prop0 = newProp('t', 'p', 1);
            d = 0.0;
            water_tp(T_K,p,d,dp,prop0);
            double h=prop0->h;
            // very very important!!!!
            prop0 = freeProp(prop0);
            return h;
        #else 
            SteamState S = freesteam_set_pT(p, T_K);
            return freesteam_h(S);
        #endif 
    }
    double cH2ONaCl::water_mu_pT(double p, double T_K)
    {
        #ifdef USE_PROST
            double d, dp, ds, dh;
            Prop *prop0;
            dp = 1.0e-8;
            ds = 1.0e-8;
            dh = 1.0e-8;
            prop0 = newProp('t', 'p', 1);
            d = 0.0;
            water_tp(T_K,p,d,dp,prop0);
            double mu=viscos(prop0);
            // very very important!!!!
            prop0 = freeProp(prop0);
            return mu;
        #else 
            SteamState S = freesteam_set_pT(p, T_K);
            return freesteam_mu(S);
        #endif 
    }
    string cH2ONaCl::checkTemperatureRange(double temperature_C)
    {
        std::string checkResult= "";
        if (temperature_C < TMIN-Kelvin || temperature_C > TMAX-Kelvin) {
            char buff[100];
            snprintf(buff, sizeof(buff), "Temperature value %.2f is out of range\n[%.1f, %.1f] ", temperature_C, TMIN - Kelvin, TMAX - Kelvin);
            checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    string cH2ONaCl::checkPressureRange(double pressure_bar)
    {
        std::string checkResult= "";
        if (pressure_bar < PMIN/1E5 || pressure_bar > PMAX/1E5) {
            char buff[100];
            snprintf(buff, sizeof(buff), "Pressure value %.2f is out of range\n[%.1f, %.1f] bar", pressure_bar, PMIN/1e5, PMAX/1e5);
            checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    string cH2ONaCl::checkSalinityRange(double salinity)
    {
        std::string checkResult= "";
        if (salinity < XMIN || salinity > XMAX) {
            char buff[100];
            snprintf(buff, sizeof(buff), "Salinity value %.4f is out of range\n[%.1f, %.1f]", salinity, XMIN, XMAX);
            checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    string cH2ONaCl::CheckRange_H(double H0, double P0, double X0)
    {
        std::string checkResult= "";

        double HMIN=1e30, HMAX=-1e30;
        double Trange[2]={TMIN, TMAX};
        double P,T,X;
        // H2ONaCl::cH2ONaCl eos;
        H2ONaCl::PROP_H2ONaCl prop;
        // for (int i = 0; i < 2; i++)
        {
        P=P0;
        for (int j = 0; j < 2; j++)
        {
            T=Trange[j];
            // for (int k = 0; k < 2; k++)
            {
            X=X0;
            prop=prop_pTX(P*1e5, T, X);
            HMIN=(prop.H<HMIN ? prop.H : HMIN);
            HMAX=(prop.H>HMAX ? prop.H : HMAX);
            }
        }
        }
        if(H0<HMIN/1000 || H0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "Enthalpy value %.2f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P = %.2f bar, X = %.4f", H0, HMIN/1000, HMAX/1000, P0, X0);
            checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    string cH2ONaCl::CheckRanges_H_PX(double HMIN0, double HMAX0, double PXrange[4])
    {
        std::string checkResult= "";
        double HMIN=1e30, HMAX=-1e30;
        double Trange[2]={TMIN, TMAX};
        double P,T,X;
        // H2ONaCl::cH2ONaCl eos;
        H2ONaCl::PROP_H2ONaCl prop;
        for (int i = 0; i < 2; i++)
        {
        P=PXrange[i];
        for (int j = 0; j < 2; j++)
        {
            T=Trange[j];
            for (int k = 0; k < 2; k++)
            {
            X=PXrange[k+2];
            prop=prop_pTX(P*1e5, T, X);
            HMIN=(prop.H<HMIN ? prop.H : HMIN);
            HMAX=(prop.H>HMAX ? prop.H : HMAX);
            }
        }
        }
        if(HMIN0<HMIN/1000 || HMIN0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "The minimum enthalpy value %.4f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P ∈ [%.2f, %.2f] bar, X ∈ [%.4f, %.4f] ",
                    HMIN0, HMIN/1000, HMAX/1000, PXrange[0], PXrange[1],PXrange[2],PXrange[3]);
            checkResult = buff;
            return checkResult;
        }
        if(HMAX0<HMIN/1000 || HMAX0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "The maximum enthalpy value %.4f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P ∈ [%.2f, %.2f] bar, X ∈ [%.4f, %.4f] ",
                    HMAX0, HMIN/1000, HMAX/1000, PXrange[0], PXrange[1],PXrange[2],PXrange[3]);checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    string cH2ONaCl::CheckRanges_H_P(double HMIN0, double HMAX0, double Prange[2], double X0)
    {
        std::string checkResult= "";
        double HMIN=1e30, HMAX=-1e30;
        double Trange[2]={TMIN, TMAX};
        double P,T,X;
        // H2ONaCl::cH2ONaCl eos;
        H2ONaCl::PROP_H2ONaCl prop;
        for (int i = 0; i < 2; i++)
        {
        P=Prange[i];
        for (int j = 0; j < 2; j++)
        {
            T=Trange[j];
            // for (int k = 0; k < 2; k++)
            {
            X=X0;
            prop=prop_pTX(P*1e5, T, X);
            HMIN=(prop.H<HMIN ? prop.H : HMIN);
            HMAX=(prop.H>HMAX ? prop.H : HMAX);
            }
        }
        }
        if(HMIN0<HMIN/1000 || HMIN0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "The minimum enthalpy value %.4f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P ∈ [%.2f, %.2f] bar, X = %.4f ",
                    HMIN0, HMIN/1000, HMAX/1000, Prange[0], Prange[1],X0);
            checkResult = buff;
            return checkResult;
        }
        if(HMAX0<HMIN/1000 || HMAX0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "The maximum enthalpy value %.4f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P ∈ [%.2f, %.2f] bar, X ∈ %.4f ",
                    HMAX0, HMIN/1000, HMAX/1000, Prange[0], Prange[1],X0);checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    string cH2ONaCl::CheckRanges_H_X(double HMIN0, double HMAX0, double Xrange[2], double P0)
    {
        std::string checkResult= "";
        double HMIN=1e30, HMAX=-1e30;
        double Trange[2]={TMIN, TMAX};
        double P,T,X;
        // H2ONaCl::cH2ONaCl eos;
        H2ONaCl::PROP_H2ONaCl prop;
        // for (int i = 0; i < 2; i++)
        {
            P=P0;
            for (int j = 0; j < 2; j++)
            {
                T=Trange[j];
                for (int k = 0; k < 2; k++)
                {
                X=Xrange[k];
                prop=prop_pTX(P*1e5, T, X);
                HMIN=(prop.H<HMIN ? prop.H : HMIN);
                HMAX=(prop.H>HMAX ? prop.H : HMAX);
                }
            }
        }
        if(HMIN0<HMIN/1000 || HMIN0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "The minimum enthalpy value %.4f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P = %.2f bar, X ∈ [%.4f, %.4f] ",
                    HMIN0, HMIN/1000, HMAX/1000, P0,Xrange[0],Xrange[1]);
            checkResult = buff;
            return checkResult;
        }
        if(HMAX0<HMIN/1000 || HMAX0>HMAX/1000)
        {
            char buff[100];
            snprintf(buff, sizeof(buff), "The maximum enthalpy value %.4f is out of range\n[%.1f, %.1f] kJ/kg corresponding to P = %.2f bar, X ∈ [%.4f, %.4f] ",
                    HMAX0, HMIN/1000, HMAX/1000, P0,Xrange[0],Xrange[1]);
            checkResult = buff;
            return checkResult;
        }
        return checkResult;
    }
    void cH2ONaCl::createTable4_Driesner2007a(TABLE4 & table4)
    {
        // Table 4 of Driesner and Heinrich(2007)
        double c[14] = {-2.36, 0.128534, -0.023707, 0.00320089, -0.000138917, 
                        1.02789E-07, -4.8376E-11, 2.36, -0.0131417, 0.00298491,
                        -0.000130114, 0, 0, -0.000488336};// c[11] and c[12] are calculated below
        double cA[11] = {1, 1.5, 2, 2.5, 3, 4, 5, 1, 2, 2.5, 3};
        // c[11] (c12 in Driesner and Heinrich(2007)) is the value of P_crit at 500 deg.C, calculated from eq. 5b.
        // c[12] (c13) is the first temperature derivative of eq. 5b at 500 deg.C
        for (size_t i = 7; i < 11; i++)
        {
            c[11] += c[i] * pow(500 - H2O::T_Critic, cA[i]); //the second term of eq. 5b in Drisner and Heinrich (2007)
            c[12] += c[i] * cA[i] * pow(500 - H2O::T_Critic, cA[i] - 1); //the first temperature derivative of eq. 5b
        }
        c[11] = H2O::P_Critic + c[11];
        double d[11] = {8E-05, 1E-05, -1.37125E-07, 9.46822E-10, -3.50549E-12, 6.57369E-15, 
                        -4.89423E-18, 7.77761E-2, 2.7042E-4, -4.244821E-07, 2.580872E-10};
        // copy 
        for(int i=0;i<14; i++)table4.c[i]    = c[i];
        for(int i=0;i<11; i++)table4.cA[i]   = cA[i];
        for(int i=0;i<11; i++)table4.d[i]    = d[i];
    }
    /**
     * - Pressure
     * 
     * \f{equation}
     *      P_{crit} = \left\{ \begin{matrix}
     *      P_{crit}^{H_2O} + \sum\limits_{n=1}^{7} c_n (T_{crit}^{H_2O} - T)^{c_nA} & , T< T_{crit}^{H_2O} \ (\text{eq. 5a})\\ \\[1ex]
     *      P_{crit}^{H_2O} + \sum\limits_{n=8}^{11} c_n (T - T_{crit}^{H_2O})^{c_nA} & , T_{crit}^{H_2O} \le T \le 500 ^{\circ}C \ (\text{eq. 5b}) \\ \\[1ex]
     *      \sum\limits_{n=12}^{14} c_n (T - 500)^{n-12} & , T > 500 ^{\circ}C \ (\text{eq. 5c})\\ \\[1ex]
     *      \end{matrix}\right.
     * \f}
     * 
     * - Salinity
     * 
     * \f{equation}
     *      X_{crit} = \left\{ \begin{matrix}
     *      \sum\limits_{i=1}^{7} d_i (T - T_{crit}^{H_2O})^{i} & , T_{crit}^{H_2O} \le T \le 600 ^{\circ}C \ (\text{eq. 7a}) \\ \\[1ex]
     *      \sum\limits_{i=8}^{11} d_i (T - 600 ^{\circ}C)^{i-8} & , 600 < T \le 1000 ^{\circ}C \ (\text{eq. 7b})\\ \\[1ex]
     *      \end{matrix}\right.
     * \f}
     * 
     * \image html Driesner_Heinrich_Fig5.png "Critical pressure." width=50%.
     * Critical pressure, figure 5 of reference \cite Driesner2007Part1. 
     * \image html Driesner_Heinrich_Fig6.png "Critical composition." width=50%.
     * Critical composition, figure 6 of reference \cite Driesner2007Part1. 
     * 
     * \image html HaliteCriticalCurves.svg "Critical pressure and composition." width=50%. 
     * Critical pressure (a,b) and composition (c,d) as function of temperature. (a,c) Full range, (b,d) the region just above the critical temperature of water
     * \warning Critical pressure calculated by Eq. 5a is different (up to 9 bar) from boiling curve of pure water. The physical meaning would be the same, but due to some numerical reason, they are not completely same.
     */
    void cH2ONaCl::P_X_Critical(double T, double& P_crit, double& X_crit)
    {
        // calculate critical pressure
        P_crit=0;
        if(T < H2O::T_Critic && T>=H2ONaCl::TMIN_C){
            for (size_t i = 0; i < 7; i++){
                P_crit += m_tab4_Driesner2007a.c[i]*pow(H2O::T_Critic - T, m_tab4_Driesner2007a.cA[i]); //eq. 5a. Note: this should give the same result as IAPWS for pure water.
            }
            P_crit+=H2O::P_Critic;
        }else if(T >= H2O::T_Critic && T <= 500){
            for (size_t i = 7; i < 11; i++){
                P_crit += m_tab4_Driesner2007a.c[i]*pow(T - H2O::T_Critic, m_tab4_Driesner2007a.cA[i]); //eq. 5b
            }
            P_crit+=H2O::P_Critic;
        }else if(T > 500 && T <= H2ONaCl::TMAX_C){
            for (size_t i = 11; i < 14; i++){
                P_crit += m_tab4_Driesner2007a.c[i]*pow(T - 500, i-11); //eq. 5c
            }
        }else
        {
            cout<<WARN_COUT<<"T: "<<T<<" out of temperature range: ["<<H2ONaCl::TMIN_C<<", "<<H2ONaCl::TMAX_C<<"]"<<endl;
        }
         // calculate critical salinity
        X_crit = 0;
        if (T>=H2O::T_Critic && T<=600){
            for (size_t i = 0; i < 7; i++){
                X_crit += m_tab4_Driesner2007a.d[i]*pow(T - H2O::T_Critic, i+1); //eq. 7a
            }
        }else if(T > 600 && T <= H2ONaCl::TMAX_C){
            for (size_t i = 7; i < 11; i++){
                X_crit += m_tab4_Driesner2007a.d[i]*pow(T - 600, i-7); //eq. 7b
            }
        }else if(T < H2ONaCl::TMIN_C || T > H2ONaCl::TMAX_C)
        {
            cout<<WARN_COUT<<"T: "<<T<<" out of temperature range: ["<<H2ONaCl::TMIN_C<<", "<<H2ONaCl::TMAX_C<<"]"<<endl;
        }
    } 
    void cH2ONaCl::P_X_Critical(std::vector<double> T, std::vector<double>& P_crit, std::vector<double>& X_crit)
    {
        P_crit.resize(T.size());
        X_crit.resize(T.size());
        for (size_t i = 0; i < T.size(); i++)
        {
            P_X_Critical(T[i], P_crit[i], X_crit[i]);
        }
        
    };
    
    void cH2ONaCl:: T_X_Critical(double P, double& T_crit, double& X_crit)
    {
        double T0 = TMAX_C, P0 = 0, P1 = 0, X_tmp = 0, dPdT=0, dT=1E-5;
        T_crit = H2O::T_Critic;
        double tol = 1E-4;
        int iter = 0;
        while (fabs(T_crit-T0)/T_crit > tol)
        {
            T0 = T_crit;
            P_X_Critical(T0, P0, X_crit); // f(x0)
            P_X_Critical(T0 + dT, P1, X_tmp);
            dPdT = (P1-P0)/dT;              //f`(x0)
            T_crit = T0 - (P0-P)/dPdT;
            // printf("%d: P=%.2f, P0=%.2f, T0=%.2f, P1=%.0f, dPdT=%.2f, Tnew = %.2f, tol=%.4f\n",iter, P, P0, T0, P1, dPdT, T, fabs(T-T0)/T);
            iter ++;
            if(iter>100)break;
        }
    };
    /**
     * The liqudius is fitted with the equation 
     * 
     * \f{equation}
     * X_{NaCl, sat}^L = \sum\limits_{i=0}^5 e_i\left(\frac{T}{T_{hm}}\right)^i
     * \f}
     * with the pressure-dependent coefficients \f$ e_i \f$ given in Table 5 of reference \cite Driesner2007Part1 and the pressure-dependent melting temperature of halite, \f$ T_{hm}\f$, calculated from NaCl::cNaCl::T_Melting (equation (1) of reference \cite Driesner2007Part1.)
     * 
     * \image html Driesner_Heinrich_Fig7.png "Liquid composition for the halite liquidus." width=50%. 
     * Liquid composition for the halite liquidus. (a) Full range temperature-pressure dependence, sets of symbols are for the same pressures as sets of lines; (b) pressure dependence at 25 \f$ ^{\circ}C \f$ (Figure 7 of reference \cite Driesner2007Part1)
     * 
     * \image html HaliteLiquidus.svg "Liquid composition for the halite liquidus calculated using H2ONaCl." width=50%. 
     * Liquid composition for the halite liquidus calculated using #H2ONaCl -> #cH2ONaCl ->#X_HaliteLiquidus. (a) Full range temperature-pressure dependence, sets of symbols are for the same pressures as sets of lines; (b) pressure dependence at 25 \f$ ^{\circ}C \f$ (Figure 7 of reference \cite Driesner2007Part1)
     *  
     */
    double cH2ONaCl::X_HaliteLiquidus(double T, double P)
    {
        const double P_squre = P*P;
        // Table 5 of Driesner and Heinrich(2007)
        double e[6] = {0.0989944 + 3.30796E-06 * P - 4.71759E-10 * P_squre, 
                       0.00947257 - 8.6646E-06 * P + 1.69417E-09 * P_squre, 
                       0.610863 - 1.51716E-05 * P + 1.1929E-08 * P_squre, 
                       -1.64994 + 0.000203441 * P - 6.46015E-08 * P_squre, 
                       3.36474 - 0.000154023 * P + 8.17048E-08 * P_squre, 
                       1
                        };
        for (size_t i = 0; i < 5; i++)
        {
            e[5] -= e[i];
        }
        const double TbyT_hm = T/m_NaCl.T_Melting(P);
        double X_Liquids = 0;
        for (size_t i = 0; i < 6; i++)
        {
            X_Liquids += e[i] * pow(TbyT_hm, i);
        }
        if(X_Liquids>1)X_Liquids=1; //ensure X_Liquids in range of [0,1]
        return X_Liquids;
    }
    std::vector<double> cH2ONaCl::X_HaliteLiquidus(std::vector<double> T, std::vector<double> P)
    {
        std::vector<double> X;
        for (size_t i = 0; i < T.size(); i++)
        {
            X.push_back(X_HaliteLiquidus(T[i],P[i]));
        }
        return X;
    };
    /**
     * The correlations for distribution coefficient \f$ K \f$ used to calculate the halite-saturated vapor composition is defined as (eq. 9 of ref. \cite Driesner2007Part1),
     * \f{equation}
     * \frac{X_{NaCl, sat}^{L, metastable}}{X_{NaCl, sat}^V} = K
     * \f}
     * A modified distribution coefficient \f$ K^{\prime} \f$ (this is actually used to compute halite-saturated vapor composition) is defined as (eq. 14 of ref. \cite Driesner2007Part1), 
     * 
     * \f{equation}
     * log_{10}K^{\prime} = log_{10}\left( \frac{x_l}{x_v/(\frac{P_{NaCl}}{P})} \right) = log_{10} \left( \frac{x_l}{x_v} \right) + log_{10} \left( \frac{P_{NaCl}}{P} \right)
     * \f}
     * where \f$ x_l \f$ is calculated from #X_HaliteLiquidus, \f$ P_{NaCl} \f$ is the boling pressure or sublimation pressure depends on temperature, it can be expressed as,
     * 
     * \f{equation}
     *      P_{NaCl} = \left\{ \begin{matrix}
     *      log_{10}(P_{NaCl, liquid}),& T > T_{triple, NaCl} \\ \\[1ex]
     *      log_{10}(P_{NaCl, halite}),& \text{else}\\ \\[1ex]
     *      \end{matrix}\right.
     * \f}
     * \f$ P_{NaCl, liquid} \f$ and \f$ P_{NaCl, halite} \f$ are calculated from NaCl::cNaCl::P_Boiling and NaCl::cNaCl::P_HaliteSublimation, respectively. And \f$ log_{10}K^{\prime} \f$ can be computed from following equation (eq. 15 of ref. \cite Driesner2007Part1),
     * 
     * \f{equation}
     * log_{10}\bar K = \frac{log_{10}K^{\prime} - log_{10}\left( X_{NaCl, sat}^L \right)_{P_{NaCl}}}{log_{10}\left( \frac{P_{NaCl}}{P_{crit}} \right) - log_{10}\left( X_{NaCl, sat}^L \right)_{P_{NaCl}}} = 1 + j_0(1-\bar P)^{j_1} + j_2(1-\bar P) + j_3(1-\bar P)^{2} - (1 + j_0 + j_2 + j_3)(1-\bar P)^{3}
     * \f}
     * where \f$ P_{crit} \f$ is calculated from #P_X_Critical, \f$ X_{NaCl, sat}^L \f$ is calculated from #X_HaliteLiquidus (T, \f$ P_{NaCl} \f$), \f$ \bar P = (P - P_{NaCl})/(P_{crit} - P_{NaCl})\f$ is the normalized pressure, \f$ j_i (i=0, 1, 2,, 3) \f$ is calculated from Table 8 of reference \cite Driesner2007Part1. 
     * 
     * \image html Driesner_Heinrich_Fig8.png "Halite-saturated vapor composition in comparison to experimental data." width=50%. 
     * Halite-saturated vapor composition in comparison to experimental data. (a) Isotherms, (b) isobars (Figure 8 of reference \cite Driesner2007Part1)
     * 
     * \image html HaliteSaturatedVaporComposition.svg "Halite-saturated vapor composition calculawithted using H2ONaCl." width=50%. 
     * Halite-saturated vapor composition calculated using #H2ONaCl -> #cH2ONaCl ->#X_VaporHaliteCoexist. 
     */
    double cH2ONaCl::X_VaporHaliteCoexist(double T, double P)
    {
        // Table 8 of Driesner and Heinrich(2007)
        double k[16] = {-0.235694, -0.188838, 0.004, 0.0552466, 0.66918, 396.848, 45, -3.2719E-07, 
                        141.699, -0.292631, -0.00139991, 1.95965E-06, -7.3653E-10, 0.904411, 0.000769766,-1.18658E-06};
        double j[4]={0,0,0,0};
        j[0] = k[0] + k[1] * exp(-k[2] * T);
        j[1] = k[4] + (k[3] - k[4]) / (1 + exp((T - k[5]) / k[6])) + k[7] * pow(T + k[8], 2.0);
        for (size_t i = 0; i < 4; i++)
        {
            j[2] += k[i+9]*pow(T,i);
        }
        for (size_t i = 0; i < 3; i++)
        {
            j[3] += k[i+13]*pow(T,i);
        }
        // cal
        double P_NaCl=0;
        if(T>NaCl::T_Triple)
        {
            P_NaCl = m_NaCl.P_Boiling(T);
        }else
        {
            P_NaCl = m_NaCl.P_Sublimation(T);
        }
        double P_crit = 0, X_crit=0;
        P_X_Critical(T,P_crit, X_crit); //calculate critic pressure
        double P_normalized = (P - P_NaCl) / (P_crit - P_NaCl); // eq. 16
        // DEBUG
        if(P_normalized>1)
        {
            // cout<<WARN_COUT<<"Normalized pressure greater than 1: "<<P_normalized<<", set it to 1"<<endl;
            P_normalized = 1;
        }
        double one_minus_P_normalized = 1 - P_normalized; //used in eq. 17
        // eq. 17
        double log10_K_overline = 1 + j[0]*pow(one_minus_P_normalized, j[1]) 
                                    + j[2]*one_minus_P_normalized 
                                    + j[3]*pow(one_minus_P_normalized, 2)
                                    - (1 + j[0] + j[2] + j[3])*pow(one_minus_P_normalized, 3);
        double log10_XL_P_NaCl = log10(X_HaliteLiquidus(T, P_NaCl)); //used in eq. 15
        
        double log10_K_prim = log10_K_overline * (log10(P_NaCl/P_crit) - log10_XL_P_NaCl) + log10_XL_P_NaCl; //eq. 15
        double log10_XLbyXV = log10_K_prim - log10(P_NaCl/P); //eq. 14
        double X_L = X_HaliteLiquidus(T, P);
        double X_V = X_L/pow(10, log10_XLbyXV); //eq. 14

        return X_V;
    }
    /**
     * \f{equation}
     * P_{VLH} = \sum\limits_{i=0}^{10}f_i\left( \frac{T}{T_{triple, NaCl}} \right)
     * \f}
     * where \f$ f_i (i=0,...,10) \f$ are computed from Table 6 of ref. \cite Driesner2007Part1, \f$ T_{triple, NaCl} \f$ is the temperature at triple point of NaCl, which is defined as #NaCl::T_Triple.
     * 
     * \image html Driesner_Heinrich_Fig9.png "Pressure at vapor + liquid + halite coexistence." width=50%. 
     * Pressure at vapor + liquid + halite coexistence. (a) Full range, (b) low temperatures, logarithmic pressure scale. (Figure 9 of reference \cite Driesner2007Part1)
     * 
     * \image html Pressure_VLH.svg "Pressure at vapor + liquid + halite coexistence calculated using H2ONaCl." width=50%. 
     * Pressure at vapor + liquid + halite coexistence calculated using #H2ONaCl -> #cH2ONaCl ->#P_VaporLiquidHaliteCoexist. 
     */
    double cH2ONaCl::P_VaporLiquidHaliteCoexist(double T)
    {
        // Table 6 of Driesner and Heinrich(2007)
        double f[11] = {0.00464, 5E-07, 16.9078, -269.148, 7632.04, -49563.6, 233119.0, -513556.0, 549708.0, -284628.0, NaCl::P_Triple};
        double P_VLH = 0;
        for (size_t i = 0; i < 10; i++) //calculate P_VLH (eq. 10) and f[10] simultaneously
        {
            f[10] -= f[i];
            P_VLH += f[i]*pow(T/NaCl::T_Triple, i);
        }
        P_VLH += f[10]*pow(T/NaCl::T_Triple, 10);
        
        return P_VLH;
    }
    void cH2ONaCl::Pmax_VaporLiquidHaliteCoexist(double& T, double& P)
    {
        // Table 6 of Driesner and Heinrich(2007)
        double f[11] = {0.00464, 5E-07, 16.9078, -269.148, 7632.04, -49563.6, 233119.0, -513556.0, 549708.0, -284628.0, NaCl::P_Triple};
        for (size_t i = 0; i < 10; i++) //calculate P_VLH (eq. 10) and f[10] simultaneously
        {
            f[10] -= f[i];
        }
        // find root
        const int degree = 9;
        Polynomial polynomial;
        std::vector<double> coefficient_vector;
        coefficient_vector.resize(degree + 1);
        double * coefficient_vector_ptr = &coefficient_vector[0];
        for (size_t i = 0; i < degree+1; i++)
        {
            coefficient_vector_ptr[i]=(i+1)*f[i+1];
        }
        polynomial.SetCoefficients(coefficient_vector_ptr, degree);

        std::vector<double> real_vector;
        std::vector<double> imag_vector;
        real_vector.resize(degree);
        imag_vector.resize(degree);
        double * real_vector_ptr = &real_vector[0];
        double * imag_vector_ptr = &imag_vector[0];
        int root_count= 0;
        if (polynomial.FindRoots(real_vector_ptr,imag_vector_ptr,&root_count) == PolynomialRootFinder::SUCCESS)
        {
            for (int i = 0; i < root_count; ++i)
            {
                if(imag_vector_ptr[i]==0)
                {
                    double root_T = real_vector_ptr[i]*NaCl::T_Triple;
                    if(root_T>H2ONaCl::TMIN_C && root_T<=NaCl::T_Triple)
                    T=root_T;
                }
            }
        }
        P = P_VaporLiquidHaliteCoexist(T);
    }
    std::vector<double> cH2ONaCl::HX_VaporLiquidHaliteCoexist(double P)
    {
        std::vector<double> HminHmaxXminXmax;
        double Hmin, Hmax, Xmin, Xmax;
        if(P>389.0 || P<H2ONaCl::PMIN)
        {
            return HminHmaxXminXmax;
        }
        // calculate T1, T2 of intersection of const P and VLH boundary surface
        vector<double> T1T2 = T_VaporLiquidHaliteCoexist(P);
        if(T1T2.size()!=2)
        {
            return HminHmaxXminXmax;
        }
        double dH=1E4;
        double Tmin = min(T1T2) -1, Tmax=max(T1T2) + 1; // prop_pTX doesn't support LVH region properties, so avoid this case at this moment, fix it later!
        double P_Pa = P*1E5;
        // 1. for Tmin point, find the exact H and X
        Xmin = Mol2Wt(X_HaliteLiquidus(Tmin, P)); // mol fraction to mass fraction
        m_prop = prop_pTX(P_Pa, Tmin+Kelvin, Xmin);
        //now the Hmin is calculated from PTX, it is close to the exact H
        Hmin = max(m_prop.H - 0.5E6, 0.1E6); // move down e.g. 0.2E6, and then move up and check phase changes
        // for this point, from low H to high H, phase changes from L or LH to LVH
        m_prop = prop_pHX(P_Pa, Hmin, Xmin);
        int iter = 0;
        while (m_prop.Region!=ThreePhase_V_L_H) //find exact Hmin in PHX space
        {
            Hmin += dH;
            m_prop = prop_pHX(P_Pa, Hmin, Xmin);
            // cout<<iter<<" Hmin: "<<Hmin<<" Xmin: "<<Xmin<<m_phaseRegion_name[m_prop.Region]<<endl;
            iter++;
            if(iter>100)return HminHmaxXminXmax;
        }
        // 2. for Tmax point, find the exact H and X
        double dX = 0.01;
        Xmax = min(Mol2Wt(X_HaliteLiquidus(Tmax, P)), 1-1E-5);
        m_prop = prop_pTX(P_Pa, Tmax+Kelvin, Xmax);
        Hmax = m_prop.H + 0.5E6; // move down e.g. 0.2E6, and then move up and check phase changes
        // for this point, from low H to high H, phase changes from L or LH to LVH
        m_prop = prop_pHX(P_Pa, Hmax, Xmax);
        iter = 0;
        while (m_prop.Region!=ThreePhase_V_L_H) //find exact Hmin in PHX space
        {
            Hmax -= dH;
            m_prop = prop_pHX(P_Pa, Hmax, Xmax);
            // cout<<iter<<" Hmax: "<<Hmax<<" Xmax: "<<Xmax<<" "<<m_phaseRegion_name[m_prop.Region]<<endl;
            iter++;
            if(iter>100)return HminHmaxXminXmax;
        }

        HminHmaxXminXmax.push_back(Hmin);
        HminHmaxXminXmax.push_back(Hmax);
        HminHmaxXminXmax.push_back(Xmin);
        HminHmaxXminXmax.push_back(Xmax);
        
        return HminHmaxXminXmax;
    }
    std::vector<double> cH2ONaCl::T_VaporLiquidHaliteCoexist(double P)
    {
        // Table 6 of Driesner and Heinrich(2007)
        // P = f0 + f1*x + f2*x^2 + ... + f10*x^{10} ==>
        // f(x)=f0 - P + f1*x + f2*x^2 + ... + f10*x^{10} = 0, where x=T/T_{triple, NaCl}
        const int degree = 10;
        double f[11] = {0.00464, 5E-07, 16.9078, -269.148, 7632.04, -49563.6, 233119.0, -513556.0, 549708.0, -284628.0, NaCl::P_Triple};
        for (size_t i = 0; i < 10; i++) 
        {
            f[10] -= f[i];
        }
        Polynomial polynomial;
        std::vector<double> coefficient_vector;
        coefficient_vector.resize(degree + 1);
        double * coefficient_vector_ptr = &coefficient_vector[0];
        for (size_t i = 0; i < degree+1; i++)
        {
            coefficient_vector_ptr[i]=f[i];
        }
        coefficient_vector_ptr[0]=coefficient_vector_ptr[0]-P;
        polynomial.SetCoefficients(coefficient_vector_ptr, degree);

        std::vector<double> real_vector;
        std::vector<double> imag_vector;
        real_vector.resize(degree);
        imag_vector.resize(degree);
        double * real_vector_ptr = &real_vector[0];
        double * imag_vector_ptr = &imag_vector[0];
        int root_count= 0;
        std::vector<double> roots_T;
        if (polynomial.FindRoots(real_vector_ptr,imag_vector_ptr,&root_count) == PolynomialRootFinder::SUCCESS)
        {
            for (int i = 0; i < root_count; ++i)
            {
                if(imag_vector_ptr[i]==0)
                {
                    double root_T = real_vector_ptr[i]*NaCl::T_Triple;
                    if(root_T>=H2ONaCl::TMIN_C && root_T<=H2ONaCl::TMAX_C)
                    roots_T.push_back(root_T);
                }
            }
        }
        
        return roots_T;
    }
    /**
     * \f{equation}
     * X_{NaCl}^{VL, liq} = X_{crit} + g_0 \sqrt{P_{crit} - P} + g_1 (P_{crit} - P) +g_2(P_{crit} - P)^2
     * \f}
     * where \f$ g_1, g_2\f$ can be found in Table 7 of ref. \cite Driesner2007Part1, \f$ P_{crit}\f$ and \f$ X_{crit}\f$ are calculated from #P_X_Critical.
     * 
     * \image html X_VaporLiquidCoexistSurface.svg "Isothermal sections of the V + L surface calculated using H2ONaCl." width=100%. 
     * Isothermal sections of the V + L surface calculated using #H2ONaCl -> #cH2ONaCl ->#X_VaporLiquidCoexistSurface_LiquidBranch for liquid branch and #X_VaporLiquidCoexistSurface_VaporBranch for vapor branch. See also Fig. 12 of ref. \cite Driesner2007Part1.
     * 
     */
    double cH2ONaCl::X_VaporLiquidCoexistSurface_LiquidBranch(double T, double P)
    {
        // Table 7 of Driesner and Heinrich(2007)
        double h[11]={0.00168486, 0.000219379, 438.58, 18.4508, -5.6765E-10, 6.73704E-06, 
                      1.44951E-07, 384.904, 7.07477, 6.06896E-05, 0.00762859};
        double g0 = 0;
        double g1 = h[1] + (h[0] - h[1]) / (1 + exp((T - h[2]) / h[3])) + h[4] * T*T;
        double g2 = h[6] + (h[5] - h[6]) / (1 + exp((T - h[7]) / h[8])) + h[9] * exp(-h[10] * T);
        double X_crit = 0, P_crit = 0;
        P_X_Critical(T, P_crit, X_crit);
        // TODO: add comment of g0 calculation formula
        double P_tmp=0, P_tmp2, X_tmp=0;
        if(T<NaCl::T_Triple)
        {
            P_tmp = P_VaporLiquidHaliteCoexist(T);
            X_tmp = X_HaliteLiquidus(T, P_tmp);
        }else
        {
            P_tmp = m_NaCl.P_Boiling(T);
            X_tmp = 1;
        }
        double X_VL_LiquidBranch = 0;
        // eq. 11
        if(T<H2O::T_Critic)
        {
            P_tmp2 = m_water.P_Boiling(T);
            g0 = (X_tmp + g1 * (P_tmp - P_tmp2) + g2 * (pow(P_crit - P_tmp2, 2.0) - pow(P_crit - P_tmp, 2.0))) / (sqrt(P_crit - P_tmp) - sqrt(P_crit - P_tmp2));
            X_VL_LiquidBranch = g0 * sqrt(P_crit - P) - g0 * sqrt(P_crit - P_tmp2) - g1 * (P_crit - P_tmp2) - g2 * pow(P_crit - P_tmp2, 2.0) + g1 * (P_crit - P) + g2 * pow(P_crit - P, 2.0);
        }else
        {
            g0 = (X_tmp - X_crit - g1 * (P_crit - P_tmp) - g2 * pow(P_crit - P_tmp, 2.0)) / sqrt(P_crit - P_tmp);
            X_VL_LiquidBranch = X_crit + g0 * sqrt(P_crit - P) + g1 * (P_crit - P) + g2 * pow(P_crit - P, 2.0);
        }

        return X_VL_LiquidBranch;
    }
    std::vector<double> cH2ONaCl::X_VaporLiquidCoexistSurface_LiquidBranch(std::vector<double> T, std::vector<double> P)
    {
        std::vector<double> X;
        for (size_t i = 0; i < T.size(); i++)
        {
            X.push_back(X_VaporLiquidCoexistSurface_LiquidBranch(T[i],P[i]));
        }
        return X;
    }
    /**
     * See also #X_VaporLiquidCoexistSurface_LiquidBranch and #X_VaporHaliteCoexist
     */
    double cH2ONaCl::X_VaporLiquidCoexistSurface_VaporBranch(double T, double P)
    {
        // Table 8 of Driesner and Heinrich(2007)
        double k[16] = {-0.235694, -0.188838, 0.004, 0.0552466, 0.66918, 396.848, 45, -3.2719E-07, 
                        141.699, -0.292631, -0.00139991, 1.95965E-06, -7.3653E-10, 0.904411, 0.000769766, -1.18658E-06};
        double j[4]={0,0,0,0};
        j[0] = k[0] + k[1] * exp(-k[2] * T);
        j[1] = k[4] + (k[3] - k[4]) / (1 + exp((T - k[5]) / k[6])) + k[7] * pow(T + k[8], 2.0);
        for (size_t i = 0; i < 4; i++)
        {
            j[2] += k[i+9]*pow(T,i);
        }
        for (size_t i = 0; i < 3; i++)
        {
            j[3] += k[i+13]*pow(T,i);
        }
        // cal
        double P_NaCl=0;
        if(T>NaCl::T_Triple)
        {
            P_NaCl = m_NaCl.P_Boiling(T);
        }else
        {
            P_NaCl = m_NaCl.P_Sublimation(T);
        }
        double X_VL_LiquidBranch = X_VaporLiquidCoexistSurface_LiquidBranch(T,P);
        double P_crit = 0, X_crit=0;
        P_X_Critical(T,P_crit, X_crit); //calculate critic pressure
        double P_normalized = (P - P_NaCl) / (P_crit - P_NaCl); // eq. 16
        // DEBUG
        if(P_normalized>1)
        {
            // cout<<WARN_COUT<<"Normalized pressure greater than 1: "<<P_normalized<<", set it to 1"<<endl;
            P_normalized = 1;
        }
        double one_minus_P_normalized = 1 - P_normalized; //used in eq. 17
        // eq. 17
        double log10_K_overline = 1 + j[0]*pow(one_minus_P_normalized, j[1]) 
                                    + j[2]*one_minus_P_normalized 
                                    + j[3]*pow(one_minus_P_normalized, 2)
                                    - (1 + j[0] + j[2] + j[3])*pow(one_minus_P_normalized, 3);

        double log10_XL_P_NaCl = log10(X_HaliteLiquidus(T, P_NaCl)); //used in eq. 15
        
        double log10_K_prim = log10_K_overline * (log10(P_NaCl/P_crit) - log10_XL_P_NaCl) + log10_XL_P_NaCl; //eq. 15

        double X_VL_VaporBranch = X_VL_LiquidBranch/pow(10, log10_K_prim) * P_NaCl/P; //eq. 14
        if(P<=P_VaporLiquidHaliteCoexist(T))
        {
            X_VL_VaporBranch = X_HaliteLiquidus(T,P)/X_VL_LiquidBranch * X_VL_VaporBranch;
        }
        return X_VL_VaporBranch;
    }
    std::vector<double> cH2ONaCl::X_VaporLiquidCoexistSurface_VaporBranch(std::vector<double> T, std::vector<double> P)
    {
        std::vector<double> X;
        for (size_t i = 0; i < T.size(); i++)
        {
            X.push_back(X_VaporLiquidCoexistSurface_VaporBranch(T[i],P[i]));
        }
        return X;
    }
    void cH2ONaCl::T_star_V_n1n2(double P, double X_NaCl, double& n1, double& n2)
    {
        double X_H2O = 1-X_NaCl;
        double P_sqrt = sqrt(P), PP = P*P, PPP=P*PP;

        double n11 = -54.2958 - 45.7623 * exp(-0.000944785 * P);
        double n1_XNaCl = 330.47 + 0.942876 * P_sqrt + 0.0817193 * P - 2.47556E-08 * PP + 3.45052E-10 * PPP; //eq. 11
        double n10 = n1_XNaCl; //eq. 10 when X_NaCl=1
        double n12 = -n11 - n10; //eq. 10 when X_NaCl=0
        n1 = n10 + n11 * X_H2O + n12 * X_H2O*X_H2O;

        double n21 = -2.6142 - 0.000239092 * P;
        double n22 = 0.0356828 + 4.37235E-06 * P + 2.0566E-09 * P*P;
        double n2_XNaCl = -0.0370751 + 0.00237723 * P_sqrt + 5.42049E-05 * P + 5.84709E-09 * PP - 5.99373E-13 * PPP; //eq. 12
        double n20 = 1 - n21 * sqrt(n22); //eq. 10 when X_NaCl=0
        double n23 = n2_XNaCl - n20 - n21 * sqrt(1 + n22); //eq. 10 when X_NaCl=1
        n2 = n20 + n21 * sqrt(X_NaCl + n22) + n23 * X_NaCl;
    }
    /**
     * \image html Driesner2007_Fig2.png "Molar volume of brine" width=25%.
     * Fig. 2 of \cite Driesner2007Part2. Graphical illustration of the principle used to derive correlations for molar volumes: the molar volume of an aqueous NaCl solution (here: 10 wt% NaCl at 1000 bar) at temperature T is identical to that of pure water at a different temperature \f$ T_V^* \f$.
     * 
     * \f{equation}
     * T_V^* = n_1 + n_2T + D(T)
     * \f}
     * where \f$ n_1, n_2\f$ are calculated from equation (9-12, 14-16) and Table 4 of \cite Driesner2007Part2. 
     * 
     * \image html V_brine_T_P1000bar.svg "Molar volume of brine calculated using swEOS." width=25%. 
     * 
     * \image html Driesner2007_Fig5.png "Examples for low temperature (a) and high temperature (b) cases" width=50%.
     * Fig. 5 of \cite Driesner2007Part2. Examples for low temperature (a) and high temperature (b) cases where the \f$ T - T_V^* \f$ correlation cannot be applied. 
     * 
     * \image html V_brine_NaCl_lowThighT.svg "Examples for low temperature (a) and high temperature (b) cases calculated by swEOS." width=50%. 
     * 
     * \bug 为什么图Fig. 5a(\cite Driesner2007Part2)不一致？因为在P=5 bar(低压)情况下，根据\f$ T_V^* \f$ 计算得到的水的密度在148 \f$^{\circ}\text{C} \f$附近存在跳跃，这是由于在此温度附近发生了相变！如下图所示。 \b 所以问题是： 为什么文献 \cite Driesner2007Part2 中的Fig. 5a能够得到一个超过150\f$^{\circ}\text{C} \f$ 的\f$ V_{sat}\f$曲线 ？
     * 
     * \image html water_rho_lowP.svg "Water density in low pressure region" width=25%.
     */
    double cH2ONaCl::T_star_V(double T, double P, double X_NaCl)
    {
        double n300 = 7606640/pow(P + 472.051, 2.0);
        double n301 = -50 - 86.1446 * exp(-0.000621128 * P);
        double n302 = 294.318 * exp(-0.00566735 * P);
        double n310 = -0.0732761 * exp(-0.0023772 * P) - 5.2948E-05 * P;
        double n311 = -47.2747 + 24.3653 * exp(-0.00125533 * P);
        double n312 = -0.278529 - 0.00081381 * P;
        double n30 = n300 * (exp(n301 * X_NaCl) - 1) + n302 * X_NaCl; //eq. 15
        double n31 = n310 * exp(n311 * X_NaCl) + n312 * X_NaCl; //eq. 16
        double D = n30 * exp(n31 * T); //eq. 14
        double n1, n2;
        T_star_V_n1n2(P,X_NaCl, n1,n2);
        return n1 + n2*T + D;
    }
    
    double cH2ONaCl::Rho_Br_for_V_extrapol(double T, double P, double X)
    {
        double T_star = T_star_V(T, P, X);
        double V_water = H2O::MolarMass / m_water.Rho(T_star, P);// m3/mol
        return (H2O::MolarMass * (1 - X) + NaCl::MolarMass * X) / V_water;
    }
    /**
     * \f{equation}
     * V_{extrapol} = o_0 + o_1 T + o_2 T^3
     * \f}
     * \todo 验证此函数的正确性，并加入brine密度计算函数
     */
    double cH2ONaCl::V_extrapol(double T, double P, double X)
    {
        //DEBUG: V_Extrapol(x_in, T_in, P_in) 
        double X_L_sat = X_HaliteLiquidus(T, P);
        double o0 = 0, o1 = 0, o2 = 0, o3=0, o4=0, o5=0, V1 = 0, V2 = 0;
        double TT=T*T;
        double TTT=T*TT;
        const double dT = 0.005;
        double V_extrapol0 = 0;
        double RoundDown_X = floor(X*1E5)/1E5, RoundDown2 = 0; //DEBUG
        double molFactor = (H2O::MolarMass * (1-X) + NaCl::MolarMass*X);
        // if (T<=200 && P<=m_water.BoilingCurve(T) && (X_L_sat - X)<0.01)
        {
            T = T_star_V(T, P, X);
            double V_L_sat = H2O::MolarMass / m_water.Rho_Liquid_Saturated(T) * 1E6; //Molar volume , cm3/mol
            double V_water = H2O::MolarMass / m_water.Rho(T, P) * 1E6;
            if (V_L_sat < V_water)
            {
                o2 = 2.0125E-07 + 3.29977E-09 * exp(-4.31279 * log10(P)) - 1.17748E-07 * log10(P) + 7.58009E-08 * pow(log10(P), 2.0); //Table 4 last row of Driesner(2007)
                V1 = H2O::MolarMass / m_water.Rho_Liquid_Saturated(T) * 1E6;//cm3/mol
                V2 = H2O::MolarMass / m_water.Rho_Liquid_Saturated(T - dT)* 1E6;//cm3/mol
                o1 = (V1 -V2)/dT - 3*o2*TT; //derivative T of eq. 17 = 0
                o0 = V1 - o1*T - o2*TTT;
                V_extrapol0 = o0 + o1*T + o2*TTT;
            }
        }
        // else if(P <= 350 && T>=600)
        // {
        //     X_L_sat = X_VaporLiquidCoexistSurface_LiquidBranch(T, P);
        //     RoundDown2 = floor(X_L_sat*1E5)/1E5; //DEBUG
        //     if(X >= X_L_sat)
        //     {
        //         double V1000 =  molFactor/ Rho_Br_for_V_extrapol(T,1000,X)* 1E6;//cm3/mol  //DEBUG: 可以简化！！！！
        //         V1 = molFactor / Rho_Br_for_V_extrapol(T,390.147,X)* 1E6;//cm3/mol
        //         V2 = molFactor / Rho_Br_for_V_extrapol(T,390.137,X)* 1E6;//cm3/mol
        //         double dVdP390 = (V1-V2)/(0.01);
        //         o4 = (V1 - V1000 + dVdP390*1609.853) / (log(1390.147 / 2000) - 2390.147 / 1390.147); //eq. 18, eq. 7 DEBUG: LogExp ????
        //         o3 = V1 - o4 * log(1390.147) - 390.147 * dVdP390 + 390.147 / 1390.147 * o4;
        //         o5 = dVdP390 - o4 /1390.147;
        //         V_extrapol0 = o3 + o4*log(P + 1000) + o5*P; //eq. 18
        //     }
        // }else
        // {
        //     V_extrapol0 = 0;
        // }
        return V_extrapol0;
    }
    double cH2ONaCl::Rho_brine(double T, double P, double X)
    {
        double T_star = T_star_V(T, P, X);
        double V_water = V_extrapol(T, P, X);
        if (V_water == 0)
        {
            V_water = H2O::MolarMass / m_water.Rho(T, P);
        }
        return (H2O::MolarMass * (1 - X) + NaCl::MolarMass * X) / V_water;
    }
    void cH2ONaCl::writeVTK_PolyLine(string filename,vector<double> X, vector<double> Y, vector<double> Z)
    {
        ofstream fpout(filename);
        if(!fpout)
        {
            cout<<"ERROR: Can not open file: "<<filename<<endl;
            exit(0);
        }
        //  write VTK head
        fpout<<"# vtk DataFile Version 2.0"<<endl;
        fpout<<"VTK 3D curve generated by swEOS(X:wt.%, T:deg.C, P:bar)"<<endl;
        fpout<<"ASCII"<<endl;
        fpout<<"DATASET POLYDATA"<<endl;
        // write points
        fpout<<"POINTS "<<X.size()<<" float"<<endl;
        for (size_t i = 0; i < X.size(); i++)
        {
            fpout<<" "<<X[i]<<" "<<Y[i]<<" "<<Z[i]<<endl;
        }
        // write vtu cells
        fpout<<"LINES 1 "<<X.size()+1<<endl;
        fpout<<X.size()<<" ";
        for (size_t i = 0; i < X.size(); i++)
        {
            fpout<<i<<" ";
        }
        fpout<<endl;
        // close file
        fpout.close();
    }
    void cH2ONaCl::writeCriticalCurve(string filename,double Tmin, double Tmax, double dT, fmtOutPutFile fmt)
    {
        // 1. Calculate 
        double P_crit=0, X_crit=0;
        vector<double> vecT, vecP, vecX;
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        for (double T = Tmin; T < Tmax; T=T+dT)
        {
            P_X_Critical(T, P_crit, X_crit);
            vecT.push_back((T-TMIN_C)/lenT);
            vecP.push_back((P_crit-PMIN)/lenP);
            vecX.push_back(Mol2Wt(X_crit)); // convert molar fraction to weight percent, [0, 1]
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_PolyLine(filename+".vtk", vecX, vecT, vecP);
            }
            break;
        
        default:
            break;
        }
    }
    void cH2ONaCl::writeNaClMeltingCurve(string filename,double Pmin, double Pmax, double dP, fmtOutPutFile fmt)
    {
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        // 1. Calculate 
        double T=0;
        vector<double> vecT, vecP, vecX;
        for (double P = Pmin; P < Pmax; P=P+dP)
        {
            T=m_NaCl.T_Melting(P);
            vecT.push_back((T-TMIN_C)/lenT);
            vecP.push_back((P-PMIN)/lenP);
            vecX.push_back(1); // weight percent NaCl, this is pure NaCl
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_PolyLine(filename+".vtk", vecX, vecT, vecP);
            }
            break;
        
        default:
            break;
        }
    }
    void cH2ONaCl::writeVaporLiquidHalite_V_L_H_Curve(string filename,double Tmin, double Tmax, fmtOutPutFile fmt, int nT)
    {
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        // 1. Calculate 
        double dT=(Tmax-Tmin)/(nT-1);
        double T, P, X_L, X_V=0, X_H=1;
        vector<double> vecT, vecP, vecX_L, vecX_V, vecX_H;
        for (size_t i = 0; i < nT; i++)
        {
            T=Tmin+i*dT;
            P = P_VaporLiquidHaliteCoexist(T);
            X_L=X_HaliteLiquidus(T,P);
            
            vecT.push_back((T-TMIN_C)/lenT);
            vecP.push_back((P-PMIN)/lenP);
            vecX_L.push_back(Mol2Wt(X_L));
            vecX_V.push_back(Mol2Wt(X_V));
            vecX_H.push_back(Mol2Wt(X_H));
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_PolyLine(filename+"_V.vtk", vecX_V, vecT, vecP);
                writeVTK_PolyLine(filename+"_L.vtk", vecX_L, vecT, vecP);
                writeVTK_PolyLine(filename+"_H.vtk", vecX_H, vecT, vecP);
            }
            break;
        
        default:
            break;
        }
    }
    void cH2ONaCl::writeH2OBoilingCurve(string filename,double Tmin, double Tmax, H2ONaCl::fmtOutPutFile fmt,int nT)
    {
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        // 1. Calculate 
        double dT=(Tmax-Tmin)/(nT-1);
        double T, P;
        vector<double> vecT, vecP, vecX;
        for (size_t i = 0; i < nT; i++)
        {
            T=Tmin+i*dT;
            P =m_water.P_Boiling(T);
            vecT.push_back((T-TMIN_C)/lenT);
            vecP.push_back((P-PMIN)/lenP);
            vecX.push_back(0);
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_PolyLine(filename+".vtk", vecX, vecT, vecP);
            }
            break;
        
        default:
            break;
        }
    }
    void cH2ONaCl::writeVTK_Triangle_Strip(string filename, vector<vector<double> > X, vector<vector<double> > Y, vector<vector<double> > Z, double scale_X, double scale_Y, double scale_Z)
    {
        ofstream fpout(filename);
        if(!fpout)
        {
            cout<<"ERROR: Can not open file: "<<filename<<endl;
            exit(0);
        }
        int np = 0;
        for (size_t i = 0; i < X.size(); i++)np+=X[i].size();
        int ncell=X.size()-1;
        //  write VTK head
        fpout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
        fpout<<"  <UnstructuredGrid>"<<endl;
        fpout<<"    <Piece NumberOfPoints=\""<<np<<"\" NumberOfCells=\""<<ncell<<"\">"<<endl;
        fpout<<"      <PointData>"<<endl;
        // fpout<<"      <PointData Scalars=\"HaliteLiquidus\">"<<endl;
        // fpout<<"        <DataArray type=\"Float64\" Name=\"HaliteLiquidus\" format=\"ascii\">"<<endl;
        // fpout<<"          ";
        //     for (size_t i = 0; i < X.size(); i++)
        //     {
        //         for (size_t j = 0; j < X[i].size(); j++)
        //         {
        //             fpout<<X[i][j]<<" ";
        //         }
        //     }
        // fpout<<"\n        </DataArray>"<<endl;
        fpout<<"      </PointData>"<<endl;
        fpout<<"      <CellData>"<<endl;
        fpout<<"      </CellData>"<<endl;
        fpout<<"      <Points>"<<endl;
        fpout<<"        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
            for (size_t i = 0; i < X.size(); i++)
            {
                for (size_t j = 0; j < X[i].size(); j++)
                {
                    fpout<<"          "<<X[i][j]*scale_X<<" "<<Y[i][j]*scale_Y<<" "<<Z[i][j]*scale_Z<<endl;
                }
            }
        fpout<<"        </DataArray>"<<endl;
        fpout<<"      </Points>"<<endl;
        fpout<<"      <Cells>"<<endl;
        fpout<<"        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
            int offset=0;
            vector<int> vec_offset;
            for (size_t i = 0; i < X.size()-1; i++)
            {
                fpout<<"          ";
                for (size_t j = 0; j < X[i].size(); j++)
                {
                    fpout<<j+offset<<" "<<j+X[i].size()+offset<<" ";
                }
                fpout<<endl;
                vec_offset.push_back(offset*2);
                offset+=X[i].size();
            }
        fpout<<"        </DataArray>"<<endl;
        fpout<<"        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;
        fpout<<"          ";
        for (size_t i = 0; i < vec_offset.size(); i++)
        {
            fpout<<vec_offset[i]<<" ";
        }
        fpout<<endl;
        fpout<<"        </DataArray>"<<endl;
        fpout<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<endl;
        fpout<<"          ";
        for (size_t i = 0; i < ncell; i++)
        {
            fpout<<6<<" "; //VTK_TRIANGLE_STRIP = 6
        }
        fpout<<endl;
        fpout<<"        </DataArray>"<<endl;
        fpout<<"      </Cells>"<<endl;
        fpout<<"    </Piece>"<<endl;
        fpout<<"  </UnstructuredGrid>"<<endl;
        fpout<<"</VTKFile>"<<endl;

        fpout.close();
    }
    void cH2ONaCl::writeVTK_Quads(string filename, vector<vector<double> > X, vector<vector<double> > Y, vector<vector<double> > Z, double scale_X, double scale_Y, double scale_Z, bool includeTwoEndsPolygon)
    {
        ofstream fpout(filename);
        if(!fpout)
        {
            cout<<"ERROR: Can not open file: "<<filename<<endl;
            exit(0);
        }
        int nPoly = X[0].size();
        int nLayers = X.size()-1;
        int nQuads = nPoly*nLayers;
        int np = X.size()*X[0].size();
        int ncell=nQuads + (includeTwoEndsPolygon==true ? 2: 0); // plus two ends or not
        //  write VTK head
        fpout<<"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<<endl;
        fpout<<"  <UnstructuredGrid>"<<endl;
        fpout<<"    <Piece NumberOfPoints=\""<<np<<"\" NumberOfCells=\""<<ncell<<"\">"<<endl;
        fpout<<"      <PointData>"<<endl;
        fpout<<"      </PointData>"<<endl;
        fpout<<"      <CellData>"<<endl;
        fpout<<"      </CellData>"<<endl;
        fpout<<"      <Points>"<<endl;
        fpout<<"        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
            for (size_t i = 0; i < X.size(); i++)
            {
                for (size_t j = 0; j < X[i].size(); j++)
                {
                    // std::setprecision(9)
                    fpout<<"          "<<X[i][j]*scale_X<<" "<<Y[i][j]*scale_Y<<" "<<Z[i][j]*scale_Z<<endl;
                }
            }
        fpout<<"        </DataArray>"<<endl;
        fpout<<"      </Points>"<<endl;
        fpout<<"      <Cells>"<<endl;
        fpout<<"        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
            vector<int> vec_offset;
            int offset = 0;
            for (size_t i = 0; i < nLayers; i++)
            {
                if (X[i].size()!=nPoly)
                {
                    printf("Error in writeVTK_Quads-> the %d layer polygon has %d nodes, but the first layer has %d nodes. This function requires node number of polygons at each layer must be the same.\n", (int)i, (int)X[i].size(), nPoly);
                    exit(0);
                }
                // each edge constructs a quad
                for (size_t j = 0; j < nPoly; j++)
                {
                    int index0 = j+i*X[i].size();
                    int index1=index0+1;
                    if(j==(nPoly-1))index1=i*X[i].size(); //the last node connect to the first node
                    fpout<<"          "<<index0<<" "<<index1<<" "<<index0+nPoly<<" "<<index1+nPoly<<endl;
                    offset += 4;
                    vec_offset.push_back(offset);
                }
            }
            if (includeTwoEndsPolygon) //include two ends polygon
            {
                for (size_t n = 0; n < 2; n++)
                {
                    fpout<<"          ";
                    for (size_t i = 0; i < nPoly; i++)
                    {
                        fpout<<i+n*nLayers*nPoly<<" ";
                    }
                    fpout<<endl;
                    offset += nPoly;
                    vec_offset.push_back(offset);
                }
            }
        fpout<<"        </DataArray>"<<endl;
        fpout<<"        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;
        fpout<<"          ";
        for (size_t i = 0; i < vec_offset.size(); i++)
        {
            fpout<<vec_offset[i]<<" ";
        }
        fpout<<endl;
        fpout<<"        </DataArray>"<<endl;
        fpout<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<endl;
        fpout<<"          ";
        // quads
        for (size_t i = 0; i < nQuads; i++)
        {
            fpout<<6<<" "; //VTK_QUAD = 9 VTK_TRIANGLE_STRIP = 6
        }
        if (includeTwoEndsPolygon)
        {
            fpout<<7<<" "<<7<<" "; //VTK_POLYGON = 7
        }
        
        fpout<<endl;
        fpout<<"        </DataArray>"<<endl;
        fpout<<"      </Cells>"<<endl;
        fpout<<"    </Piece>"<<endl;
        fpout<<"  </UnstructuredGrid>"<<endl;
        fpout<<"</VTKFile>"<<endl;
        fpout.close();

        // // write gmsh geo file
        // string filename_gmsh=filename+".geo";
        // FILE* fp_gmsh=fopen(filename_gmsh.c_str(),"w");
        // fprintf(fp_gmsh,"SetFactory(\"OpenCASCADE\");\n");
        // fprintf(fp_gmsh,"lc=1;\n");
        // int index_point = 0, index_surface=0;
        // for (size_t i = 0; i < X.size(); i++)
        // {
        //     fprintf(fp_gmsh,"Point(%d) = {%f, %f, %f, lc};\n",index_point, X[i][0]*scale_X, Y[i][0]*scale_Y, Z[i][0]*scale_Z);
        //     index_point++;
        //     for (size_t j = 1; j < X[i].size(); j++)
        //     {
        //         fprintf(fp_gmsh,"Point(%d) = {%f, %f, %f, lc};\n",index_point, X[i][j]*scale_X, Y[i][j]*scale_Y, Z[i][j]*scale_Z);
        //         fprintf(fp_gmsh,"Line(%d) = {%d, %d};\n",index_point,index_point-1, index_point);

        //         index_point++;
        //     }
        //     fprintf(fp_gmsh,"Line(%d) = {%d, %zu};\n",index_point,index_point-1,index_point-X[i].size());
        //     fprintf(fp_gmsh,"Curve Loop(%zu) = {%zu:%d};\n",i,index_point-X[i].size()+1,index_point);
        //     index_surface++;
        // }
        // fprintf(fp_gmsh,"vol=newv;ThruSections(vol) = {%d:%d};\n",0,index_surface-1);
        // fclose(fp_gmsh);
    }
    void cH2ONaCl::writeHaliteLiquidusSurface(string filename,double Tmin, double Tmax, double Pmax, fmtOutPutFile fmt, int nT, int nP)
    {
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        double T=0,P=0, X=0, Pmin=0;
        double dT=(Tmax-Tmin)/(nT-1);
        double T_NaClMelting = 0;
        vector<vector<double> > vec2T, vec2P, vec2X;
        // 1. calculate
        for (size_t i = 0; i < nT; i++)
        {
            vector<double> vecT, vecP, vecX;
            T=Tmin+i*dT;
            Pmin = P_VaporLiquidHaliteCoexist(T);
            double dP=(Pmax-Pmin)/(nP-1);
            if (T<NaCl::T_Triple)
            {
                for (size_t j = 0; j < nP; j++)
                {
                    P = Pmin + j*dP;
                    X=X_HaliteLiquidus(T,P);
                    vecT.push_back((T-TMIN_C)/lenT);
                    vecP.push_back((P-PMIN)/lenP);
                    vecX.push_back(Mol2Wt(X)); // convert molar fraction to weight percent, [0, 1]
                }
            }else
            {
                for (size_t j = 0; j < nP; j++)
                {
                    P = Pmin + j*dP;
                    T = m_NaCl.T_Melting(P);
                    X=X_HaliteLiquidus(T,P);
                    vecT.push_back((T-TMIN_C)/lenT);
                    vecP.push_back((P-PMIN)/lenP);
                    vecX.push_back(Mol2Wt(X)); // convert molar fraction to weight percent, [0, 1]
                }
                
            }
            vec2T.push_back(vecT);
            vec2P.push_back(vecP);
            vec2X.push_back(vecX);
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_Triangle_Strip(filename+".vtu", vec2X, vec2T, vec2P);
            }
            break;

        default:
            break;
        }
    }
    void cH2ONaCl::writeVaporLiquidHaliteCoexistSurface(string filename, double Tmin, double Tmax, double dT, fmtOutPutFile fmt)
    {
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        vector<vector<double> > vec2T, vec2P, vec2X;
        int nX=101;
        double Xmin=0, Xmax=1;
        double dx=(Xmax - Xmin)/(nX-1);
        // 1. calculate
        double P=0;
        for (double T = Tmin; T <=Tmax; T=T+dT)
        {
            vector<double> vecT, vecP, vecX;
            P = P_VaporLiquidHaliteCoexist(T);
            for (size_t i = 0; i < nX; i++)
            {
                double X=Xmin+i*dx;
                vecT.push_back((T-TMIN_C)/lenT);
                vecP.push_back((P-PMIN)/lenP);
                vecX.push_back(Mol2Wt(X)); // convert molar fraction to weight percent, [0, 1]
            }
            vec2T.push_back(vecT);
            vec2P.push_back(vecP);
            vec2X.push_back(vecX);
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_Triangle_Strip(filename+".vtu", vec2X, vec2T, vec2P);
            }
            break;

        default:
            break;
        }
    }
    void cH2ONaCl::writeVaporLiquidCoexistSurface(string filename, double Tmin, double Tmax, H2ONaCl::fmtOutPutFile fmt,int nT, int nP)
    {
        const double lenT=TMAX_C-TMIN_C, lenP=PMAX-PMIN, lenX=XMAX-XMIN;
        vector<vector<double> > vec2T, vec2P, vec2T_v, vec2P_v, vec2X_l, vec2X_v;
        double dT=(Tmax-Tmin)/(nT-1);
        double T, P, X, Pmin, Pmax, X_crit, X_l, X_v;
        // 1. calculate
        for (size_t i = 0; i < nT; i++)
        {
            T = Tmin + i*dT;
            Pmin = PMIN;
            if(T<NaCl::T_Triple)Pmin = P_VaporLiquidHaliteCoexist(T);
            if(T>H2O::T_Critic)
            {
                P_X_Critical(T,Pmax, X_crit);
            }else
            {
                Pmax = m_water.P_Boiling(T);
            }
            double dP = (Pmax-Pmin)/(nP-1);
            vector<double> vecT, vecP, vecX_l, vecX_v;
            int n_dP_refine = 2;
            for (size_t j = 0; j < nP-n_dP_refine; j++)
            {
                P = Pmin + j*dP;
                X_l = X_VaporLiquidCoexistSurface_LiquidBranch(T,P);
                X_v = X_VaporLiquidCoexistSurface_VaporBranch(T,P);
                vecT.push_back((T-TMIN_C)/lenT);
                vecP.push_back((P-PMIN)/lenP);
                vecX_l.push_back(Mol2Wt(X_l)); // convert molar fraction to weight percent, [0, 1]
                vecX_v.push_back(Mol2Wt(X_v)); // convert molar fraction to weight percent, [0, 1]
            }
            // refine where close to top
            int np_refine = nP;
            double Pmin_refine=Pmin + (nP-n_dP_refine-1)*dP;
            double dP_refine=(Pmax-Pmin_refine)/(np_refine-1);
            for (size_t j = 0; j < np_refine; j++)
            {
                P = Pmin_refine + j*dP_refine;
                X_l = X_VaporLiquidCoexistSurface_LiquidBranch(T,P);
                X_v = X_VaporLiquidCoexistSurface_VaporBranch(T,P);
                vecT.push_back((T-TMIN_C)/lenT);
                vecP.push_back((P-PMIN)/lenP);
                vecX_l.push_back(Mol2Wt(X_l)); // convert molar fraction to weight percent, [0, 1]
                vecX_v.push_back(Mol2Wt(X_v)); // convert molar fraction to weight percent, [0, 1]
            }
            
            vec2T.push_back(vecT);
            vec2P.push_back(vecP);
            vec2X_l.push_back(vecX_l);
            // if temperature too low, vapor branch area is too small.
            if(T>50)
            {
                vec2T_v.push_back(vecT);
                vec2P_v.push_back(vecP);
                vec2X_v.push_back(vecX_v);
            }
        }
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_Triangle_Strip(filename+"_V.vtu", vec2X_v, vec2T_v, vec2P_v);
                writeVTK_Triangle_Strip(filename+"_L.vtu", vec2X_l, vec2T, vec2P);
            }
            break;

        default:
            break;
        }
    }
    PhaseRegion cH2ONaCl::findPhaseRegion(const double T_c, const double P_bar, const double X_wt, double& Xl_all, double& Xv_all)
    {
        return findRegion(T_c, P_bar*1E5, Xwt2Xmol(X_wt), Xl_all, Xv_all);
    }
    PhaseRegion cH2ONaCl::findPhaseRegion(const double T_c, const double P_bar, const double X_wt)
    {
        double Xl_all, Xv_all;
        return findRegion(T_c, P_bar*1E5, Xwt2Xmol(X_wt), Xl_all, Xv_all);
    }
    PhaseRegion cH2ONaCl::findPhaseRegion_pTX(double p_Pa, double T_K, double X_wt)
    {
        double Xl_all, Xv_all;
        return findRegion(T_K - 273.15, p_Pa, Xwt2Xmol(X_wt), Xl_all, Xv_all);
    }
    void cH2ONaCl::writePhaseSurface_XHP(double scale_X, double scale_H, double scale_P, string outpath, H2ONaCl::fmtOutPutFile fmt, int nP)
    {
        double Pmin = PMIN; // bar
        double Pmax=1000; //PMAX;
        
        double dP = (Pmax - Pmin)/(nP-1);
        int nT=100;
        // for case of P < Pmax_VLH 
        std::vector<std::vector<double> > vec2P, vec2H_part1, vec2X_part1, vec2H_part2, vec2X_part2, vec2X_VH, vec2H_VH, vec2P_VH;
        std::vector<std::vector<double> > vec2P_LH_part1, vec2P_LH_part2, vec2H_LH_part1, vec2X_LH_part1, vec2H_LH_part2, vec2X_LH_part2;
        std::vector<std::vector<double> > vec2P_VL_part1, vec2H_VL_part1, vec2X_VL_part1, vec2P_VL_part2, vec2H_VL_part2, vec2X_VL_part2;
        std::vector<std::vector<double> > vec2P_L_part1, vec2H_L_part1, vec2X_L_part1;
        std::vector<std::vector<double> > vec2P_L_part2, vec2H_L_part2, vec2X_L_part2;
        // for case of P > Pmax_VLH 
        std::vector<std::vector<double> > vec2P_LH_highP, vec2H_LH_highP, vec2X_LH_highP;
        std::vector<std::vector<double> > vec2P_VL_highP, vec2H_VL_highP, vec2X_VL_highP;
        std::vector<std::vector<double> > vec2P_L, vec2H_L, vec2X_L;
        double P0, P0_Pa;
        for (size_t i = 0; i < nP; i++)
        {
            P0 = Pmin+dP*i;
            P0_Pa = P0*1E5;
            std::vector<double> HminHmaxXminXmax = HX_VaporLiquidHaliteCoexist(P0);
            if (HminHmaxXminXmax.size()==4)
            {
                std::vector<double> vecP(3, P0), vecX_part1(3), vecX_part2(3), vecH_part1(3), vecH_part2(3), vecX_VH(4), vecP_VH(4, P0), vecH_VH(4);
                std::vector<double> vecP_LH_part1(nT+2,P0), vecH_LH_part1(nT+2), vecX_LH_part1(nT+2);
                std::vector<double> vecP_LH_part2(nT+1,P0), vecH_LH_part2(nT+1), vecX_LH_part2(nT+1);
                // for P=P0, the VLH region in salinity-enthalpy space is a triangle
                // part 1
                H2ONaCl::PROP_H2ONaCl prop1 = prop_pHX(P0_Pa, HminHmaxXminXmax[0], HminHmaxXminXmax[2]);
                vecH_part1[0] = prop1.H_v;
                vecH_part1[1] = prop1.H_l;
                vecH_part1[2] = prop1.H_h;
                vecX_part1[0] = 0;
                vecX_part1[1] = HminHmaxXminXmax[2];
                vecX_part1[2] = 1.0;
                // V+H
                vecX_VH[0] = 0; vecX_VH[1] = 1;
                vecH_VH[0] = prop1.H_v; vecH_VH[1] = prop1.H_h;
                // L+H
                double dT = (prop1.T -1 - TMIN_C)/(nT-1);
                for (size_t k = 0; k < nT; k++)
                {
                    double T_LH = TMIN_C + k*dT;
                    double X_LH=Mol2Wt(X_HaliteLiquidus(T_LH, P0));
                    m_prop = prop_pTX(P0_Pa, T_LH+Kelvin, X_LH);
                    vecH_LH_part1[k]=m_prop.H;
                    vecX_LH_part1[k]=X_LH;
                }
                vecH_LH_part1[nT]=prop1.H_h; vecX_LH_part1[nT]=1.0; //corner 1
                m_prop = prop_pTX(P0_Pa, TMIN_K, 1.0);
                vecH_LH_part1[nT+1]=m_prop.H; vecX_LH_part1[nT+1]=1.0; //corner 2
                vec2X_LH_part1.push_back(vecX_LH_part1);
                vec2H_LH_part1.push_back(vecH_LH_part1);
                vec2P_LH_part1.push_back(vecP_LH_part1);
                // VL
                double T_crit = 0, X_crit = 0;
                T_X_Critical(P0, T_crit, X_crit); //Critical point
                T_crit+=1E-5; //avoid critical point
                // -------- refine close to critical point -----------
                int nT_refine = (int)(nT/3.0*2);
                int nT_normal = nT - nT_refine;
                double Tmax_VL = prop1.T - 1;
                double deltaT_refine = (Tmax_VL - T_crit)*0.03;
                double T_crit_bigger = T_crit + deltaT_refine; 
                double dT_refine = deltaT_refine/(nT_refine - 1);
                double dT_normal = (Tmax_VL - T_crit_bigger)/(nT_normal - 1);
                std::vector<double> vecT(nT);
                for (size_t k = 0; k < nT_refine; k++)
                {
                    vecT[k] = T_crit + k*dT_refine;
                }
                for (size_t k = 0; k < nT_normal; k++)
                {
                    vecT[k+nT_refine] = T_crit_bigger + k*dT_normal;
                }
                // --------------
                std::vector<double> vecP_VL_part1(nT*2,P0), vecH_VL_part1(nT*2), vecX_VL_part1(nT*2);
                dT = (prop1.T -1 - T_crit)/(nT - 1); //Avoid too closing to VLH point
                H2ONaCl::PROP_H2ONaCl prop_VL_L, prop_VL_V;
                double X_VL_L, X_VL_V, T_VL_L, T_VL_V;
                for (size_t k = 0; k < nT; k++)
                {
                    T_VL_L = vecT[k]; //T_crit + dT*k;
                    T_VL_V = vecT[k];
                    // Liquid branch
                    X_VL_L = Mol2Wt(X_VaporLiquidCoexistSurface_LiquidBranch(T_VL_L,P0));
                    X_VL_L = ((X_VL_L<=0 || isnan(X_VL_L)) ? 1E-8 : X_VL_L);
                    prop_VL_L = prop_pTX(P0_Pa, T_VL_L + Kelvin, X_VL_L);
                    vecH_VL_part1[nT + k] = prop_VL_L.H;
                    vecX_VL_part1[nT + k] = X_VL_L;
                    // printf("T=%.2f P=%.2f X=%.2f H=%.2f\n", T_VL, P0, X_VL_L, prop_VL_L.H);
                    // vapor branch
                    X_VL_V = Mol2Wt(X_VaporLiquidCoexistSurface_VaporBranch(T_VL_V,P0));
                    X_VL_V = ((X_VL_V<=0 || isnan(X_VL_V)) ? 1E-8 : X_VL_V);
                    prop_VL_V = prop_pTX(P0_Pa, T_VL_V + Kelvin, X_VL_V);
                    vecH_VL_part1[nT - k - 1] = prop_VL_V.H;
                    vecX_VL_part1[nT - k - 1] = X_VL_V;
                }
                vec2P_VL_part1.push_back(vecP_VL_part1);
                vec2H_VL_part1.push_back(vecH_VL_part1);
                vec2X_VL_part1.push_back(vecX_VL_part1);
                // Single phase liquid
                std::vector<double> vecH_L_part1, vecX_L_part1;
                for (size_t k = 0; k < vecH_VL_part1.size(); k++)
                {
                    vecH_L_part1.push_back(vecH_VL_part1[k]);
                    vecX_L_part1.push_back(vecX_VL_part1[k]);
                }
                for (int k = vecH_LH_part1.size()-3; k > -1; k--)
                {
                    vecH_L_part1.push_back(vecH_LH_part1[k]);
                    vecX_L_part1.push_back(vecX_LH_part1[k]);
                }
                m_prop = prop_pTX(P0_Pa, TMIN_K, XMIN+1E-5);
                vecH_L_part1.push_back(m_prop.H);
                vecX_L_part1.push_back(0);

                vector<double> vecP_L_part1(vecH_L_part1.size(), P0);
                vec2H_L_part1.push_back(vecH_L_part1);
                vec2X_L_part1.push_back(vecX_L_part1);
                vec2P_L_part1.push_back(vecP_L_part1);

                // part 2
                H2ONaCl::PROP_H2ONaCl prop2 = prop_pHX(P0_Pa, HminHmaxXminXmax[1], HminHmaxXminXmax[3]);
                vecH_part2[0] = prop2.H_v;
                vecH_part2[1] = prop2.H_l;
                vecH_part2[2] = prop2.H_h;
                vecX_part2[0] = 0;
                vecX_part2[1] = HminHmaxXminXmax[3];
                vecX_part2[2] = 1.0;
                // V+H region
                vecX_VH[2] = 1; vecX_VH[3] = 0;
                vecH_VH[2] = prop2.H_h; vecH_VH[3] = prop2.H_v;

                // L+H
                double T_VLH2 = prop2.T;
                dT = (NaCl::T_Triple - T_VLH2)/(nT-1);
                for (size_t k = 0; k < nT; k++)
                {
                    double T_LH = T_VLH2 +1 + k*dT;
                    double X_LH=Mol2Wt(X_HaliteLiquidus(T_LH, P0));
                    m_prop = prop_pTX(P0_Pa, T_LH+Kelvin, X_LH);
                    // printf("T=%.2f P=%.2f X=%.2f H=%.2f, T2=%.2f\n", T_LH, P0, X_LH, m_prop.H,T_VLH2);
                    vecH_LH_part2[k]=m_prop.H;
                    vecX_LH_part2[k]=X_LH;
                }
                vecH_LH_part2[nT]=prop2.H_h; vecX_LH_part2[nT]=1.0; //corner 1
                vec2X_LH_part2.push_back(vecX_LH_part2);
                vec2H_LH_part2.push_back(vecH_LH_part2);
                vec2P_LH_part2.push_back(vecP_LH_part2);
                // V+L
                // calculate maximum H and corresponding X according to maximum T
                H2ONaCl::PROP_H2ONaCl prop_XminTmax = prop_pTX(P0_Pa, TMAX_K, XMIN+1E-8);
                std::vector<double> vecP_VL_part2(nT+2,P0), vecH_VL_part2(nT+2), vecX_VL_part2(nT+2);
                dT = (TMAX_C - prop2.T)/(nT - 1); 
                vecH_VL_part2[0] = prop2.H_v;
                vecX_VL_part2[0] = 0;
                for (size_t k = 0; k < nT; k++)
                {
                    T_VL_L = prop2.T + dT*k;
                    // Liquid branch
                    X_VL_L = Mol2Wt(X_VaporLiquidCoexistSurface_LiquidBranch(T_VL_L,P0));
                    // X_VL_L = ((X_VL_L<=0 || isnan(X_VL_L)) ? 1E-8 : X_VL_L);
                    prop_VL_L = prop_pTX(P0_Pa, T_VL_L + Kelvin, X_VL_L);
                    vecH_VL_part2[k+1] = prop_VL_L.H;
                    vecX_VL_part2[k+1] = X_VL_L;
                }
                vecH_VL_part2[nT+1] = prop_XminTmax.H;
                vecX_VL_part2[nT+1] = 0;
                vec2P_VL_part2.push_back(vecP_VL_part2);
                vec2H_VL_part2.push_back(vecH_VL_part2);
                vec2X_VL_part2.push_back(vecX_VL_part2);
                // Single phase liquid
                std::vector<double> vecH_L_part2, vecX_L_part2;
                for (size_t k = 1; k < vecH_VL_part2.size()-1; k++)
                {
                    vecH_L_part2.push_back(vecH_VL_part2[k]);
                    vecX_L_part2.push_back(vecX_VL_part2[k]);
                }
                m_prop = prop_pTX(P0_Pa, TMAX_K, XMAX-1E-5);
                vecH_L_part2.push_back(m_prop.H);
                vecX_L_part2.push_back(1);
                for (int k = (vecH_LH_part2.size()-2); k > -1; k--)
                {
                    vecH_L_part2.push_back(vecH_LH_part2[k]);
                    vecX_L_part2.push_back(vecX_LH_part2[k]);
                }
                vector<double> vecP_L_part2(vecH_L_part2.size(), P0);
                vec2H_L_part2.push_back(vecH_L_part2);
                vec2X_L_part2.push_back(vecX_L_part2);
                vec2P_L_part2.push_back(vecP_L_part2);

                // =======
                vec2P.push_back(vecP);

                vec2X_part1.push_back(vecX_part1);
                vec2H_part1.push_back(vecH_part1);

                vec2X_part2.push_back(vecX_part2);
                vec2H_part2.push_back(vecH_part2);

                vec2P_VH.push_back(vecP_VH);
                vec2X_VH.push_back(vecX_VH);
                vec2H_VH.push_back(vecH_VH);
            }
            else
            {
                // L+H
                vector<double> vecH_LH(nT+1), vecX_LH(nT+1), vecP_LH(nT+1, P0);
                double dT = (NaCl::T_Triple - TMIN_C)/(nT-1);
                for (size_t k = 0; k < nT; k++)
                {
                    double T_LH = TMIN_C + k*dT;
                    double X_LH=Mol2Wt(X_HaliteLiquidus(T_LH, P0));
                    m_prop = prop_pTX(P0_Pa, T_LH+Kelvin, X_LH);
                    vecH_LH[k]=m_prop.H;
                    vecX_LH[k]=X_LH;
                }
                m_prop = prop_pTX(P0_Pa, TMIN_K, XMAX);
                vecH_LH[nT]=m_prop.H;   vecX_LH[nT]=XMAX;
                vec2H_LH_highP.push_back(vecH_LH);
                vec2X_LH_highP.push_back(vecX_LH);
                vec2P_LH_highP.push_back(vecP_LH);
                // VL
                double T_crit = 0, X_crit = 0;
                T_X_Critical(P0, T_crit, X_crit); //Critical point
                T_crit+=1E-5;
                std::vector<double> vecP_VL_highP(nT*2,P0), vecH_VL_highP(nT*2), vecX_VL_highP(nT*2);
                int nT_refine = (int)(nT/3.0*2);
                int nT_normal = nT - nT_refine;
                double deltaT_refine = (TMAX_C - T_crit)*0.1;
                double T_crit_bigger = T_crit + deltaT_refine; 
                double dT_refine = deltaT_refine/(nT_refine - 1);
                double dT_normal = (TMAX_C - T_crit_bigger)/(nT_normal - 1);
                std::vector<double> vecT(nT);
                for (size_t k = 0; k < nT_refine; k++)
                {
                    vecT[k] = T_crit + k*dT_refine;
                }
                for (size_t k = 0; k < nT_normal; k++)
                {
                    vecT[k+nT_refine] = T_crit_bigger + k*dT_normal;
                }
                H2ONaCl::PROP_H2ONaCl prop_VL_L, prop_VL_V;
                double X_VL_L, X_VL_V, T_VL;
                for (size_t k = 0; k < nT; k++)
                {
                    T_VL = vecT[k];
                    // Liquid branch
                    X_VL_L = Mol2Wt(X_VaporLiquidCoexistSurface_LiquidBranch(T_VL,P0));
                    prop_VL_L = prop_pTX(P0_Pa, T_VL + Kelvin, X_VL_L);
                    vecH_VL_highP[nT + k] = prop_VL_L.H;
                    vecX_VL_highP[nT + k] = X_VL_L;
                    // vapor branch
                    X_VL_V = Mol2Wt(X_VaporLiquidCoexistSurface_VaporBranch(T_VL,P0));
                    X_VL_V = ((X_VL_V<=0 || isnan(X_VL_V)) ? 1E-8 : X_VL_V);
                    prop_VL_V = prop_pTX(P0_Pa, T_VL + Kelvin, X_VL_V);
                    vecH_VL_highP[nT - k -1] = prop_VL_V.H;
                    vecX_VL_highP[nT - k -1] = X_VL_V;
                }
                vec2P_VL_highP.push_back(vecP_VL_highP);
                vec2H_VL_highP.push_back(vecH_VL_highP);
                vec2X_VL_highP.push_back(vecX_VL_highP);
                // Single phase liquid
                std::vector<double> vecH_L, vecX_L;
                for (size_t k = 0; k < vecP_VL_highP.size(); k++)
                {
                    vecH_L.push_back(vecH_VL_highP[k]);
                    vecX_L.push_back(vecX_VL_highP[k]);
                }
                m_prop = prop_pTX(P0_Pa, TMAX_K, XMAX-1E-5);
                vecH_L.push_back(m_prop.H);
                vecX_L.push_back(1);
                for (int k = vecH_LH.size()-2; k > -1; k--)
                {
                    vecH_L.push_back(vecH_LH[k]);
                    vecX_L.push_back(vecX_LH[k]);
                }
                m_prop = prop_pTX(P0_Pa, TMIN_K, XMIN+1E-5);
                vecH_L.push_back(m_prop.H);
                vecX_L.push_back(0);

                vector<double> vecP_L(vecH_L.size(), P0);
                vec2H_L.push_back(vecH_L);
                vec2X_L.push_back(vecX_L);
                vec2P_L.push_back(vecP_L);
            }
            
        }
        cout<<"write to vtk files"<<endl;
        // 2. Write
        switch (fmt)
        {
        case fmt_vtk:
            {
                writeVTK_Quads(outpath+"/VLH_lowH.vtu", vec2X_part1, vec2H_part1, vec2P, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/VLH_highH.vtu", vec2X_part2, vec2H_part2, vec2P, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/VH.vtu", vec2X_VH, vec2H_VH, vec2P_VH, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/LH_lowH.vtu", vec2X_LH_part1, vec2H_LH_part1, vec2P_LH_part1, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/LH_highH.vtu", vec2X_LH_part2, vec2H_LH_part2, vec2P_LH_part2, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/VL_lowH.vtu", vec2X_VL_part1, vec2H_VL_part1, vec2P_VL_part1, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/VL_highH.vtu", vec2X_VL_part2, vec2H_VL_part2, vec2P_VL_part2, scale_X, scale_H, scale_P);
                // for case of P>Pmax_VLH
                writeVTK_Quads(outpath+"/LH_highP.vtu", vec2X_LH_highP, vec2H_LH_highP, vec2P_LH_highP, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/VL_highP.vtu", vec2X_VL_highP, vec2H_VL_highP, vec2P_VL_highP, scale_X, scale_H, scale_P);
                writeVTK_Quads(outpath+"/L_highP.vtu", vec2X_L, vec2H_L, vec2P_L, scale_X, scale_H, scale_P, false);
                writeVTK_Quads(outpath+"/L_lowH.vtu", vec2X_L_part1, vec2H_L_part1, vec2P_L_part1, scale_X, scale_H, scale_P, false);
                writeVTK_Quads(outpath+"/L_highH.vtu", vec2X_L_part2, vec2H_L_part2, vec2P_L_part2, scale_X, scale_H, scale_P, false);
            }
            break;

        default:
            break;
        }
    }

    void cH2ONaCl::createLUT_2D_TPX(std::string type, double xy_min[2], double xy_max[2], double constZ, int min_level, int max_level, string filename_vtu)
    {
        clock_t start = clock();
        STATUS("Creating 2D lookup table ...");
        const int dim =2;
        m_lut_PTX_2D = new LookUpTableForest_2D (xy_min, xy_max, constZ, LOOKUPTABLE_FOREST::CONST_X, LOOKUPTABLE_FOREST::EOS_SPACE_TPX, max_level, this);
        // refine
        m_lut_PTX_2D->set_min_level(min_level);
        m_lut_PTX_2D->refine(refine_uniform);
        // parallel refine
        #pragma omp parallel //shared(n)
        {
            #pragma omp single
            {
                printf("Do refinement using %d threads.\n", m_num_threads);
                m_lut_PTX_2D->refine(RefineFunc_PTX_consX);
            }
        }
        // m_lut_PTX_2D->refine(RefineFunc_PTX_consX);
        STATUS_time("Lookup table refinement done", (clock() - start)/m_num_threads);
        
    }

    void cH2ONaCl::save_lut_to_vtk(string filename)
    {
        
        m_lut_PTX_2D->write_to_vtk(filename);
    }

    void cH2ONaCl::save_lut_to_binary(string filename)
    {
        m_lut_PTX_2D->write_to_binary(filename);
    }

    void cH2ONaCl::destroyLUT()
    {
        if(m_lut_PTX_2D) 
        {
            m_lut_PTX_2D->destory();
            delete m_lut_PTX_2D;
            m_lut_PTX_2D = NULL;
        }
        if(m_lut_PTX_3D) 
        {
            m_lut_PTX_3D->destory();
            delete m_lut_PTX_3D;
            m_lut_PTX_3D = NULL;
        }
    }
    void cH2ONaCl::createLUT_2D_TPX(std::string type, double xmin, double xmax, double ymin, double ymax, double constZ, int min_level, int max_level, string filename_vtu)
    {
        double xy_min[2] = {xmin, ymin};
        double xy_max[2] = {xmax, ymax};
        createLUT_2D_TPX(type, xy_min, xy_max, constZ, min_level,max_level, filename_vtu);
    }

    template<int dim>
    void cH2ONaCl::interp_quad_prop(LOOKUPTABLE_FOREST::Quadrant<dim,LOOKUPTABLE_FOREST::FIELD_DATA<dim> > *targetLeaf, H2ONaCl::PROP_H2ONaCl& prop, const double xyz[dim])
    {
        double physical_length[dim]; //physical length of the quad
        double coeff[dim][2];
        double values_at_vertices[m_lut_PTX_2D->m_num_children];
        m_lut_PTX_2D->get_quadrant_physical_length(targetLeaf->level, physical_length);
        get_coeff_bilinear<dim> (targetLeaf->xyz, physical_length, xyz, coeff);
        // Rho
        for (int i = 0; i < m_lut_PTX_2D->m_num_children; i++){
            values_at_vertices[i] = targetLeaf->user_data->prop_point[i].Rho;
        }
        bilinear_cal<dim>(coeff, values_at_vertices, prop.Rho);
        // H
        for (int i = 0; i < m_lut_PTX_2D->m_num_children; i++){
            values_at_vertices[i] = targetLeaf->user_data->prop_point[i].H;
        }
        bilinear_cal<dim>(coeff, values_at_vertices, prop.H);

        // phase region
        prop.Region = targetLeaf->user_data->phaseRegion_cell;
    }

    LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > * cH2ONaCl::searchLUT_2D_PTX(H2ONaCl::PROP_H2ONaCl& prop, double x, double y)
    {
        LOOKUPTABLE_FOREST::Quadrant<2,LOOKUPTABLE_FOREST::FIELD_DATA<2> > *targetLeaf = NULL;
        m_lut_PTX_2D->searchQuadrant(targetLeaf, x, y, m_lut_PTX_2D->m_constZ);
        if(targetLeaf->user_data->need_refine)
        {
            prop = prop_pTX(y, x, m_lut_PTX_2D->m_constZ);
        }
        else
        {
            double xy[2] = {x, y};
            interp_quad_prop(targetLeaf, prop, xy);
        }
        return targetLeaf;
    }

    // for python API
    H2ONaCl::PROP_H2ONaCl cH2ONaCl::searchLUT_2D_PTX(double x, double y)
    {
        H2ONaCl::PROP_H2ONaCl prop;
        searchLUT_2D_PTX(prop, x, y);
        return prop;
    }

    void cH2ONaCl::loadLUT_PTX(string filename)
    {
        int dim = LOOKUPTABLE_FOREST::get_dim_from_binary(filename);
        switch (dim)
        {
        case 2:
            m_lut_PTX_2D = (LookUpTableForest_2D*)(new LookUpTableForest_2D(filename, this));
            break;
        case 3:
            m_lut_PTX_3D = (LookUpTableForest_3D*)(new LookUpTableForest_3D(filename, this));
            break;
        default:
            break;
        }
        
    }
}
