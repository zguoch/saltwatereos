#include "H2O.H"

namespace H2O
{
    cH2O::cH2O(/* args */)
    {
        //load table 6.2 data in \cite wagner2002iapws
        LoadTable62(m_Table62);
        m_isHighAccuracy = true; //using high accuracy scheme in Rho function
    }
    
    cH2O::~cH2O()
    {
    }
    /**
     * @brief Load Table 6.2 of \cite wagner2002iapws.
     * 
     * @param m_Table62 
     */
    void cH2O::LoadTable62(Table62& m_Table62)
    {
        // initialize as zero
        for (size_t i = 0; i < 56; i++)
        {
            m_Table62.c[i] = 0;
            m_Table62.d[i] = 0;
            m_Table62.t[i] = 0;
            m_Table62.n[i] = 0;
            m_Table62.alpha[i] = 0;
            m_Table62.beta[i] = 0;
            m_Table62.gamma[i] = 0;
            m_Table62.epsilon[i] = 0;
            m_Table62.a[i] = 0;
            m_Table62.b[i] = 0;
            m_Table62.A[i] = 0;
            m_Table62.B[i] = 0;
            m_Table62.C[i] = 0;
            m_Table62.D[i] = 0;
        }
        // c
        for (size_t i = 7; i < 22; i++)m_Table62.c[i]=1;
        for (size_t i = 22; i < 42; i++)m_Table62.c[i]=2;
        for (size_t i = 42; i < 46; i++)m_Table62.c[i]=3;
        m_Table62.c[46]=4;
        for (size_t i = 47; i < 51; i++)m_Table62.c[i]=6;
        
        // d
        m_Table62.d[0]=1;     m_Table62.d[1]=1;     m_Table62.d[2]=1;     m_Table62.d[3]=2;     m_Table62.d[4]=2;     m_Table62.d[5]=3;     m_Table62.d[6]=4;
        m_Table62.d[7]=1;     m_Table62.d[8]=1;     m_Table62.d[9]=1;     m_Table62.d[10]=2;    m_Table62.d[11]=2;    m_Table62.d[12]=3;    m_Table62.d[13]=4;
        m_Table62.d[14]=4;    m_Table62.d[15]=5;    m_Table62.d[16]=7;    m_Table62.d[17]=9;    m_Table62.d[18]=10;   m_Table62.d[19]=11;   m_Table62.d[20]=13;
        m_Table62.d[21]=15;   m_Table62.d[22]=1;    m_Table62.d[23]=2;    m_Table62.d[24]=2;    m_Table62.d[25]=2;    m_Table62.d[26]=3;    m_Table62.d[27]=4;
        m_Table62.d[28]=4;    m_Table62.d[29]=4;    m_Table62.d[30]=5;    m_Table62.d[31]=6;    m_Table62.d[32]=6;    m_Table62.d[33]=7;    m_Table62.d[34]=9;
        m_Table62.d[35]=9;    m_Table62.d[36]=9;    m_Table62.d[37]=9;    m_Table62.d[38]=9;    m_Table62.d[39]=10;   m_Table62.d[40]=10;   m_Table62.d[41]=12;
        m_Table62.d[42]=3;    m_Table62.d[43]=4;    m_Table62.d[44]=4;    m_Table62.d[45]=5;    m_Table62.d[46]=14;   m_Table62.d[47]=3;    m_Table62.d[48]=6;
        m_Table62.d[49]=6;    m_Table62.d[50]=6;    m_Table62.d[51]=3;    m_Table62.d[52]=3;    m_Table62.d[53]=3;

        // t
        m_Table62.t[0]= -0.5; m_Table62.t[1]= 0.875; m_Table62.t[2]= 1;   m_Table62.t[3]= 0.5; m_Table62.t[4]= 0.75; m_Table62.t[5]= 0.375; 
        m_Table62.t[6]= 1;    m_Table62.t[7]= 4;     m_Table62.t[8]= 6;   m_Table62.t[9]= 12; m_Table62.t[10]= 1;   m_Table62.t[11]= 5; 
        m_Table62.t[12]= 4;   m_Table62.t[13]= 2;    m_Table62.t[14]= 13; m_Table62.t[15]= 9;  m_Table62.t[16]= 3;   m_Table62.t[17]= 4; 
        m_Table62.t[18]= 11;  m_Table62.t[19]= 4;    m_Table62.t[20]= 13; m_Table62.t[21]= 1;  m_Table62.t[22]= 7;   m_Table62.t[23]= 1; 
        m_Table62.t[24]= 9;   m_Table62.t[25]= 10;   m_Table62.t[26]= 10; m_Table62.t[27]= 3;  m_Table62.t[28]= 7;   m_Table62.t[29]= 10; 
        m_Table62.t[30]= 10;  m_Table62.t[31]= 6;    m_Table62.t[32]= 10; m_Table62.t[33]= 10; m_Table62.t[34]= 1;   m_Table62.t[35]= 2; 
        m_Table62.t[36]= 3;   m_Table62.t[37]= 4;    m_Table62.t[38]= 8;  m_Table62.t[39]= 6;  m_Table62.t[40]= 9;   m_Table62.t[41]= 8; 
        m_Table62.t[42]= 16;  m_Table62.t[43]= 22;   m_Table62.t[44]= 23; m_Table62.t[45]= 23; m_Table62.t[46]= 10;  m_Table62.t[47]= 50; 
        m_Table62.t[48]= 44;  m_Table62.t[49]= 46;   m_Table62.t[50]= 50; m_Table62.t[51]= 0;  m_Table62.t[52]= 1;   m_Table62.t[53]= 4; 

        // n
        m_Table62.n[0]= 0.012533547935523;      m_Table62.n[1]= 7.8957634722828; 
        m_Table62.n[2]= -8.7803203303561;       m_Table62.n[3]= 0.31802509345418; 
        m_Table62.n[4]= -0.26145533859358;      m_Table62.n[5]= -0.0078199751687981; 
        m_Table62.n[6]= 0.0088089493102134;     m_Table62.n[7]= -0.66856572307965; 
        m_Table62.n[8]= 0.20433810950965;       m_Table62.n[9]= -6.6212605039687E-05; 
        m_Table62.n[10]= -0.19232721156002;     m_Table62.n[11]= -0.25709043003438; 
        m_Table62.n[12]= 0.16074868486251;      m_Table62.n[13]= -0.040092828925807; 
        m_Table62.n[14]= 3.9343422603254E-07;   m_Table62.n[15]= -7.5941377088144E-06; 
        m_Table62.n[16]= 0.00056250979351888;   m_Table62.n[17]= -1.5608652257135E-05; 
        m_Table62.n[18]= 1.1537996422951E-09;   m_Table62.n[19]= 3.6582165144204E-07; 
        m_Table62.n[20]= -1.3251180074668E-12; m_Table62.n[21]= -6.2639586912454E-10; 
        m_Table62.n[22]= -0.10793600908932;     m_Table62.n[23]= 0.017611491008752; 
        m_Table62.n[24]= 0.22132295167546;      m_Table62.n[25]= -0.40247669763528; 
        m_Table62.n[26]= 0.58083399985759;      m_Table62.n[27]= 0.0049969146990806; 
        m_Table62.n[28]= -0.031358700712549;    m_Table62.n[29]= -0.74315929710341; 
        m_Table62.n[30]= 0.4780732991548;       m_Table62.n[31]= 0.020527940895948; 
        m_Table62.n[32]= -0.13636435110343;     m_Table62.n[33]= 0.014180634400617; 
        m_Table62.n[34]= 0.0083326504880713;    m_Table62.n[35]= -0.029052336009585; 
        m_Table62.n[36]= 0.038615085574206;     m_Table62.n[37]= -0.020393486513704; 
        m_Table62.n[38]= -0.0016554050063734;   m_Table62.n[39]= 0.0019955571979541; 
        m_Table62.n[40]= 0.00015870308324157;   m_Table62.n[41]= -1.638856834253E-05; 
        m_Table62.n[42]= 0.043613615723811;     m_Table62.n[43]= 0.034994005463765; 
        m_Table62.n[44]= -0.076788197844621;    m_Table62.n[45]= 0.022446277332006; 
        m_Table62.n[46]= -6.2689710414685E-05;  m_Table62.n[47]= -5.5711118565645E-10; 
        m_Table62.n[48]= -0.19905718354408;     m_Table62.n[49]= 0.31777497330738; 
        m_Table62.n[50]= -0.11841182425981;     m_Table62.n[51]= -31.306260323435; 
        m_Table62.n[52]= 31.546140237781;       m_Table62.n[53]= -2521.3154341695; 
        m_Table62.n[54]= -0.14874640856724;     m_Table62.n[55]= 0.31806110878444; 
        // alpha
        m_Table62.alpha[51]= 20; 
        m_Table62.alpha[52] = 20; 
        m_Table62.alpha[53] = 20;
        // beta
        m_Table62.beta[51] = 150;
        m_Table62.beta[52] = 150;
        m_Table62.beta[53] = 250;
        m_Table62.beta[54] = 0.3;
        m_Table62.beta[55] = 0.3;
        // gamma
        m_Table62.gamma[51] = 1.21;
        m_Table62.gamma[52] = 1.21;
        m_Table62.gamma[53] = 1.25;
        // epislon
        m_Table62.epsilon[51] = 1;
        m_Table62.epsilon[52] = 1;
        m_Table62.epsilon[53] = 1;
        // a
        m_Table62.a[54] = 3.5;
        m_Table62.a[55] = 3.5;
        // b
        m_Table62.b[54] = 0.85;
        m_Table62.b[55] = 0.95;
        // A
        m_Table62.A[54] = 0.32;
        m_Table62.A[55] = 0.32;
        // B
        m_Table62.B[54] = 0.2;
        m_Table62.B[55] = 0.2;
        // C
        m_Table62.C[54] = 28;
        m_Table62.C[55] = 32;
        // D
        m_Table62.D[54] = 700;
        m_Table62.D[55] = 800;

    }
    /**
     * \f{equation}
     * ln\left( \frac{P_{boil}}{P_c} \right) = \frac{T_c + 273.15}{T + 273.15}(a_1\theta + a_2\theta^{1.5} + a_3\theta^{3} + a_4\theta^{3.5} + a_5\theta^{4} + a_6\theta^{7.5}), 
     * \f} 
     * where \f$ \theta = 1- T/T_c \f$, \f$ T_c\f$ and \f$ P_c\f$ are the critical temperature and pressure, respectively. They are defined as #T_Critic and #P_Critic respectively. Coefficients \f$ a_i (i=1, ..., 6) \f$ can be found in page 399 of reference \cite wagner2002iapws.
     */
    double cH2O::P_Boiling(double T)
    {
        double T_K =T;
        if(T_K==0)T_K=0.01;
        T_K = T_K + Kelvin;
        double T_inv = 1 - T_K / T_Critic_K;
        double a[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
        double P_boil = 0;

        P_boil = exp((T_Critic + Kelvin) / T_K * (a[0] * T_inv + a[1] * pow(T_inv,1.5) + a[2] * pow(T_inv, 3.0) + a[3] * pow(T_inv, 3.5) + a[4] * pow(T_inv,4.0) + a[5] * pow(T_inv ,7.5))) * P_Critic;
        return P_boil;
    }
    double cH2O::T_Boiling(double P)
    {
        if(P>P_Critic)return NAN;

        double T0 = TMIN, P0 = 0, P1 = 0, dPdT=0, dT=1E-5;
        double T_boiling = (T_Critic - TMIN)/2.0;
        double tol = 1E-4;
        int iter = 0;
        while (fabs(T_boiling-T0)/T_boiling > tol)
        {
            T0 = T_boiling;
            P0=P_Boiling(T0); // f(x0)
            P1=P_Boiling(T0 + dT);
            dPdT = (P1-P0)/dT;              //f`(x0)
            T_boiling = T0 - (P0-P)/dPdT;
            if(T_boiling>T_Critic)T_boiling=T_Critic-dT*10;
            // printf("%d: P=%.2f, P0=%.2f, T0=%.2f, P1=%.0f, dPdT=%.2f, Tnew = %.2f, tol=%.4f\n",iter, P, P0, T0, P1, dPdT, T_boiling, fabs(T_boiling-T0)/T_boiling);
            iter ++;
            if(iter>100)break;
        }
        return T_boiling;
    }
    /**
     * \f{equation}
     * \frac{\rho_{Liquid, sat}}{\rho_c} = 1 + b_1\theta^{1/3} + b_2\theta^{2/3} + b_3\theta^{5/3} + b_4\theta^{16/3} + b_5\theta^{43/3} + b_6\theta^{110/3}
     * \f}
     * where \f$ \theta = 1- T/T_c \f$, \f$ \rho_c\f$ is the critical density defined as #Rho_Critic. Coefficients \f$ b_i (i=1, ..., 6) \f$ can be found in page 399 of reference \cite wagner2002iapws.
     */
    double cH2O::Rho_Liquid_Saturated(double T)
    {
        double T_K = T;
        if(T_K==0)T_K=0.01;
        T_K = T_K + Kelvin;
        double T_inv = 1 - T_K / T_Critic_K;
        double b[6] = {1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45};
        double rho_L_sat =0;
        rho_L_sat = (1 
                    + b[0] * pow(T_inv, 1/3.0) 
                    + b[1] * pow(T_inv, 2/3.0) 
                    + b[2] * pow(T_inv, 5/3.0) 
                    + b[3] * pow(T_inv, 16/3.0) 
                    + b[4] * pow(T_inv, 43/3.0) 
                    + b[5] * pow(T_inv, 110/3.0)) * Rho_Critic;
        return rho_L_sat;
    }
    /**
     * \f{equation}
     * \frac{\rho_{Vapor, sat}}{\rho_c} = c_1\theta^{2/6} + c_2\theta^{4/6} + c_3\theta^{8/6} + c_4\theta^{18/6} + c_5\theta^{37/6} + c_6\theta^{71/6}
     * \f}
     * where \f$ \theta = 1- T/T_c \f$, \f$ \rho_c\f$ is the critical density defined as #Rho_Critic. Coefficients \f$ c_i (i=1, ..., 6) \f$ can be found in page 399 of reference \cite wagner2002iapws.
     */
    double cH2O::Rho_Vapor_Saturated(double T)
    {
        double T_K = T;
        if(T_K==0)T_K=0.01;
        T_K = T_K + Kelvin;
        double T_inv = 1 - T_K / T_Critic_K;
        double c[6] = {-2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581, -63.9201063};
        double rho_V_sat =0;
        rho_V_sat = exp(c[0] * pow(T_inv, 1/3.0) + 
                        c[1] * pow(T_inv, 2/3.0) + 
                        c[2] * pow(T_inv, 4/3.0) + 
                        c[3] * pow(T_inv, 3.0) + 
                        c[4] * pow(T_inv, 37/6.0) + 
                        c[5] * pow(T_inv, 71/6.0)) * Rho_Critic;
        return rho_V_sat;
    }
    /**
     * \f{equation}
     \left\{ \begin{matrix}
     *      p = \rho^2 \left( \frac{\partial f}{\partial \rho} \right)_T\\ \\[1ex]
     *      \frac{p(\delta, \tau)}{\rho R T} = 1 + \delta \phi_{\delta}^{r}\\ \\[1ex]
     *      \end{matrix}\right.
     * \f}
     * where \f$ f = f(\rho, T) \f$ is the <a href="https://en.wikipedia.org/wiki/Helmholtz_free_energy">specific Helmholtz free energy</a>, \f$ \phi(\delta, \tau) = f(\rho, T)/(RT) = \phi^o(\delta, \tau) + \phi^r(\delta, \tau) \f$,  \f$ \phi^o \f$ and \f$ \phi^r \f$ represent ideal-gas part and residual part, respectively. 
     * \f$ \delta = \rho/\rho_c, \tau = T_c/T \f$ with \f$ \rho_c, T_c, R\f$ being the critical density, critical temperature and specific gas constant, and defined as #Rho_Critic, #T_Critic and #R_const, respectively. **Note that** unit of temperautre in the equations is [\f$ K \f$], pressure is [Pa].
     */
    double cH2O::Pressure_T_Rho(double T, double Rho)
    {
        double T_K = T + Kelvin;
        double delta = Rho/Rho_Critic;
        double tau = T_Critic_K/T_K;
        double pressure = Rho*R_const*T_K*(1 + delta*Phi_r_delta(delta, tau))/100.0; //note that R_const unit is kJ/kg/K, and return pressure with unit of bar.
        
        return pressure;
    }
    /**
     * \sa #SublimationCurve and #MeltingCurve
     */
    double cH2O::BoilingCurve(double T)
    {
        double Rho_Vapor_sat = Rho_Vapor_Saturated(T);
        return Pressure_T_Rho(T, Rho_Vapor_sat);
    }
    void cH2O::BoilingCurve(std::vector<double> T, std::vector<double>& p_boiling)
    {
        p_boiling.resize(T.size());
        double Rho_Vapor_sat = 0;
        for (size_t i = 0; i < T.size(); i++)
        {
            Rho_Vapor_sat = Rho_Vapor_Saturated(T[i]);
            p_boiling[i] = Pressure_T_Rho(T[i], Rho_Vapor_sat);
        }
    }
    double cH2O::Rho(double T, double P)
    {//DEBUG: Rho_Water# in water_prop.vb
        double T_K = T + Kelvin;
        double Rho1 = 0, Rho2 = 0;
        if(T_K <= T_Critic_K)
        {
            // double RoundDown_P = floor(P*1000)/1000.0;//DEBUG
            // double RoundDown_BoilingCurve = floor(BoilingCurve(T)*1000)/1000.0;//DEBUG
            if(P <= BoilingCurve(T))
            {
                Rho1 = 1E-6;
                Rho2 = Rho_Vapor_Saturated(T) + 1;
            }else
            {
                Rho1 = Rho_Liquid_Saturated(T) -1;
                Rho2 = 1701; 
            }
        }
        else
        {
            Rho1 = 1E-6;
            Rho2 = 1701;
        }
        // Bisscectional aproximation first
        const int iterMax = 1000;
        double Tol = 1; //initial tolerance
        double Rho_approx = 0, P_Rho_approx = 0, P_Rho1 = 0;
        int n=0;
        while (n<=iterMax)
        {
            Rho_approx = (Rho1 + Rho2)/2.0;
            if(fabs(Rho2 - Rho1)/2.0 < Tol)
            {
                Rho_approx = Rho2;
                n = iterMax;
            }
            n++;
            P_Rho_approx = Pressure_T_Rho(T, Rho_approx) - P; //bar
            P_Rho1 = Pressure_T_Rho(T, Rho1)- P; //bar
            if ((P_Rho_approx * P_Rho1)>0)
            {
                Rho1 = Rho_approx;
            }else
            {
                Rho2 = Rho_approx;
            }
        }
        if(Rho1 < 0)Rho1 = 0;
        if(Rho2 < 0)Rho2 = 0;
        //Newthon method to compare desired P_in with pressure, and to get accuracy defined as either below 0.01 or 1e-10 bar
        if(m_isHighAccuracy)Tol = 1E-5;
        else Tol = 1E-4;
        
        Rho1 = Rho_approx;
        Rho2 = Rho1 - Tol;
        n=0;
        bool endLoop = false;
        double P_Rho2 = 0, dPdRho = 0;
        while (!endLoop)
        {
            Rho_approx = Rho2;
            P_Rho1 = Pressure_T_Rho(T, Rho1) - P; //bar
            P_Rho2 = Pressure_T_Rho(T, Rho2) - P; //bar
            dPdRho = (P_Rho1 - P_Rho2)/(Rho1 - Rho2);
            Rho1 = Rho1 - P_Rho1/dPdRho;
            P_Rho1 = P_Rho1 + P;
            Rho2 = Rho1 - Tol;
            n++;
            if (n>=10000)//double check, to avoid dead loop!
            {
                Rho1 = NAN; 
                endLoop = true;
            }
            if (fabs(1 - P/P_Rho1)<=1E-8 || fabs(Rho_approx - Rho1)<Tol)
            {
                endLoop = true;
            }
        }
        return Rho1;
    }
    /**
     * \image html Wagner2002_Fig2.1.png "The phase boundary curves of water in a p-T diagram." width=25%.
     * Fig. 2.1 of \cite wagner2002iapws. The phase–boundary curves of water in a p–T diagram. The sublimation curve psubl and the several melting curves pm plotted in bold correspond to the lower limit of the range of validity of IAPWS-95 with regard to temperature. The dashed line corresponds to the range of validity of IAPWS-95 with regard to pressure.
     * 
     * \image html water_PhaseDiagram.svg "Phase diagram calculated by swEOS" width=50%.
     * where Sublimation curve, boiling curve and melting curve are calculated by #SublimationCurve, #BoilingCurve and #MeltingCurve, respectively. The triple point and critical point are defined as (#T_Critic, #P_Critic) and (#T_Triple, #P_Triple), respectively.
     */
    double cH2O::SublimationCurve(double T)
    {
        double T_K = T + Kelvin;
        double Pn = 0.000611657; //MPa
        double Tn = 273.16; // K
        double theta = T_K / Tn;
        return 10 * Pn * exp(-13.928169 * (1 - pow(theta, -1.5)) + 34.7078238 * (1 - pow(theta, -1.25))); //MPa -> bar
    }
    /**
     * The melting curve is separated into 5 segments, they are ice I, III, V, VI, VII. The corresponding equation can be written as(see equation 2.16-2.20 of \cite wagner2002iapws),
     * \f{equation}
     * \left\{ \begin{matrix}
     *  p_{m,ice\ I} = ..., &T\in [251.165, T_{triple}] K \\ \\[1ex]
     *  p_{m,ice\ III} = ..., &T\in [251.165, 256.164] K \\ \\[1ex]
     *  p_{m,ice\ V} = ..., &T\in [256.164, 273.31] K \\ \\[1ex]
     *  p_{m,ice\ VI} = ..., &T\in [273.31, 355] K \\ \\[1ex]
     *  p_{m,ice\ VII} = ..., &T\in [355, 715] K \\ \\[1ex]
     * \end{matrix}\right.
     * \f}
     * \image html Wagner2002_Fig2.1.png "The phase boundary curves of water in a p-T diagram." width=25%.
     * Fig. 2.1 of \cite wagner2002iapws. The phase–boundary curves of water in a p–T diagram. The sublimation curve psubl and the several melting curves pm plotted in bold correspond to the lower limit of the range of validity of IAPWS-95 with regard to temperature. The dashed line corresponds to the range of validity of IAPWS-95 with regard to pressure.
     * 
     * \sa #SublimationCurve and #BoilingCurve
     */
    double cH2O::MeltingCurve(double T, bool isIceI)
    {
        double T_K = T + Kelvin;
        double P_m =0;
        if(isIceI && T_K>=T_K_ice_min[iceI] && T_K<T_K_ice_max[iceI]) //segment I
        {
            double Tn = 273.15;
            double Pn = 0.000611657;
            double theta = T_K / Tn;
            P_m = Pn * (1 - 0.626 * 1E6 * (1 - pow(theta, -3)) + 0.197135 * 1E6 * (1 - pow(theta, 21.2))); //MPa
        }else if(T_K>=T_K_ice_min[iceIII] && T_K<T_K_ice_max[iceIII]) //ice III
        {
            double Tn = 251.165;
            double Pn = 209.9;
            double theta = T_K / Tn;
            P_m = Pn * (1 - 0.0295252 * (1 - pow(theta, 60.0)));
        }else if(T_K>=T_K_ice_min[iceV] && T_K<T_K_ice_max[iceV]) //ice V
        {
            double Tn = 256.164;
            double Pn = 350.1;
            double theta = T_K / Tn;
            P_m = Pn * (1 - 1.18721 * (1 - pow(theta, 8.0)));
        }else if(T_K>=T_K_ice_min[iceVI] && T_K<T_K_ice_max[iceVI]) //ice VI
        {
            double Tn = 273.31;
            double Pn = 632.4;
            double theta = T_K / Tn;
            P_m = Pn * (1 - 1.07476 * (1 - pow(theta, 4.6)));
        }else if(T_K>=T_K_ice_min[iceVII] && T_K<T_K_ice_max[iceVII]) //ice VII
        {
            double Tn = 355;
            double Pn = 2216;
            double theta = T_K / Tn;
            P_m = Pn * exp(1.73683 * (1 - 1.0/theta) - 0.0544606 * (1 - pow(theta, 5.0)) + 0.806106 * 1E-07 * (1 - pow(theta, 22.0)));
        }else
        {
            double Tn = 273.15;
            double Pn = 0.000611657;
            double theta = T_K / Tn;
            P_m = Pn * (1 - 0.626 * 1E6 * (1 - pow(theta, -3)) + 0.197135 * 1E6 * (1 - pow(theta, 21.2))); //MPa
        }
        
        
        return P_m * 10.0; //MPa to bar
    }
    double cH2O::Phi_r_delta(double delta, double tau)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
        for (size_t i = 0; i < 7; i++)
        {
            sum1 += m_Table62.n[i] * m_Table62.d[i] * pow(delta, m_Table62.d[i] -1) * pow(tau, m_Table62.t[i]);
        }
        for (size_t i = 7; i < 51; i++)
        {
            sum2 += m_Table62.n[i] * exp(-pow(delta,m_Table62.c[i])) * 
                    (pow(delta, m_Table62.d[i] -1) * pow(tau, m_Table62.t[i]) * 
                    (m_Table62.d[i] - m_Table62.c[i]*pow(delta, m_Table62.c[i])));
        }
        for (size_t i = 51; i < 54; i++)
        {
            sum3 += m_Table62.n[i] * pow(delta,m_Table62.d[i]) * pow(tau, m_Table62.t[i]) * exp(-m_Table62.alpha[i]*pow(delta - m_Table62.epsilon[i], 2.0) - m_Table62.beta[i]*pow(tau - m_Table62.gamma[i], 2.0) ) * (m_Table62.d[i]/delta - 2*m_Table62.alpha[i]*(delta - m_Table62.epsilon[i]) );
        }
        double psi = 0, theta = 0, Delta = 0, dpsiddelta = 0, dDeltaddelta = 0, dDeltabiddelta = 0, delta_minus_one_squre = 0, tau_minus_one_squre=0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);
            tau_minus_one_squre = pow(tau - 1, 2.0);
            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*tau_minus_one_squre);

            theta = (1-tau) + m_Table62.A[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i]);

            Delta = pow(theta, 2.0) + m_Table62.B[i]*pow(delta_minus_one_squre, m_Table62.a[i]);
            dpsiddelta = -2*m_Table62.C[i]*(delta - 1) * psi;

            dDeltaddelta = 2*(delta - 1)*(theta * m_Table62.A[i]/m_Table62.beta[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i] - 1) + m_Table62.B[i]*m_Table62.a[i]*pow(delta_minus_one_squre, m_Table62.a[i]-1));
            if(Delta == 0)
            {
                dDeltabiddelta = 0;
            }else
            {
                dDeltabiddelta = dDeltaddelta*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-1);
            }
            
            sum4 += m_Table62.n[i]*(pow(Delta, m_Table62.b[i]) * (psi + delta*dpsiddelta) + dDeltabiddelta*delta*psi);
        }
        
        return sum1+sum2+sum3+sum4;
    }
    double cH2O::Phi_r_deltadelta(double delta, double tau)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
        for (size_t i = 0; i < 7; i++)
        {
            sum1 += m_Table62.n[i] * m_Table62.d[i] * (m_Table62.d[i] - 1) * pow(delta,m_Table62.d[i] - 2) * pow(tau,m_Table62.t[i]);
        }
        for (size_t i = 7; i < 51; i++)
        {
            sum2 += m_Table62.n[i] * exp(-pow(delta,m_Table62.c[i])) * (pow(delta,m_Table62.d[i]-2) * pow(tau,m_Table62.t[i]) * ((m_Table62.d[i] - m_Table62.c[i] * pow(delta,m_Table62.c[i])) * (m_Table62.d[i] - 1 - m_Table62.c[i] * pow(delta,m_Table62.c[i])) - pow(m_Table62.c[i],2.0) * pow(delta,m_Table62.c[i])));
        }
        for (size_t i = 51; i < 54; i++)
        {
            sum3 += m_Table62.n[i] * pow(tau,m_Table62.t[i]) * exp(-m_Table62.alpha[i] * pow(delta - m_Table62.epsilon[i], 2.0) - m_Table62.beta[i] * pow(tau - m_Table62.gamma[i], 2.0)) * (-2 * m_Table62.alpha[i] * pow(delta,m_Table62.d[i]) + 4 * pow(m_Table62.alpha[i], 2.0) * pow(delta,m_Table62.d[i]) * pow(delta-m_Table62.epsilon[i],2.0) - 4 * m_Table62.d[i] * m_Table62.alpha[i] * pow(delta,m_Table62.d[i]-1) * (delta - m_Table62.epsilon[i]) + m_Table62.d[i] * (m_Table62.d[i] - 1) * pow(delta, m_Table62.d[i]-2));
        }
        double psi = 0, theta = 0, Delta = 0, dpsiddelta = 0, d2Psiddeltadelta=0, dDeltaddelta = 0, d2Deltaddeltadelta=0, dDeltabiddelta = 0,d2DeltaBIddeltadelta=0, delta_minus_one_squre = 0, tau_minus_one_squre=0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);
            tau_minus_one_squre = pow(tau - 1, 2.0);
            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*tau_minus_one_squre);
            theta = (1-tau) + m_Table62.A[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i]);
            Delta = pow(theta, 2.0) + m_Table62.B[i]*pow(delta_minus_one_squre, m_Table62.a[i]);

            dDeltaddelta = 2*(delta - 1)*(theta * m_Table62.A[i]/m_Table62.beta[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i] - 1) + m_Table62.B[i]*m_Table62.a[i]*pow(delta_minus_one_squre, m_Table62.a[i]-1));

            dDeltabiddelta = dDeltaddelta*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-1);

            dpsiddelta = -2*m_Table62.C[i]*(delta - 1) * psi;

            d2Psiddeltadelta = (2 * m_Table62.C[i] * pow(delta - 1, 2.0) - 1) * 2 * m_Table62.C[i] * psi;
            
            d2Deltaddeltadelta = 1 / (delta - 1) * dDeltaddelta + pow(delta - 1, 2) * (4 * m_Table62.B[i] * m_Table62.a[i] * (m_Table62.a[i] - 1) * pow(pow(delta - 1, 2), m_Table62.a[i]-2) + 2 * pow(m_Table62.A[i], 2)* pow(m_Table62.beta[i], -2) * pow(pow(pow(delta - 1, 2), 0.5/m_Table62.beta[i]-1), 2) + m_Table62.A[i] * theta * 4 / m_Table62.beta[i] * (0.5/ m_Table62.beta[i] - 1) * pow(pow(delta - 1, 2), 0.5/m_Table62.beta[i]-2));
            
            d2DeltaBIddeltadelta = m_Table62.b[i] * (pow(Delta, m_Table62.b[i]-1) * d2Deltaddeltadelta + (m_Table62.b[i] - 1) * pow(Delta, m_Table62.b[i]-2) * pow(dDeltaddelta,2));

            sum4 += m_Table62.n[i]*(pow(Delta,m_Table62.b[i]) * (2 * dpsiddelta + delta * d2Psiddeltadelta) + 2 * dDeltabiddelta * (psi + delta * dpsiddelta) + d2DeltaBIddeltadelta * psi * delta);
        }
        
        return sum1+sum2+sum3+sum4;
    }
    double cH2O::Phi_r_deltatau(double delta, double tau)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
        for (size_t i = 0; i < 7; i++)
        {
            sum1 += m_Table62.n[i] * m_Table62.d[i] * m_Table62.t[i] * pow(delta, m_Table62.d[i]-1.0) * pow(tau, m_Table62.t[i]-1.0);
        }
        for (size_t i = 7; i < 51; i++)
        {
            sum2 += m_Table62.n[i] * m_Table62.t[i] * pow(delta, m_Table62.d[i]-1.0) * pow(tau, m_Table62.t[i]-1.0) * (m_Table62.d[i] - m_Table62.c[i] * pow(delta, m_Table62.c[i])) * exp(-pow(delta,m_Table62.c[i]));
        }
        for (size_t i = 51; i < 54; i++)
        {
            sum3 += m_Table62.n[i] * pow(delta, m_Table62.d[i]) * pow(tau,m_Table62.t[i]) * exp(-m_Table62.alpha[i] * pow(delta - m_Table62.epsilon[i], 2.0) - m_Table62.beta[i] * pow(tau - m_Table62.gamma[i], 2.0)) * (m_Table62.d[i] / delta - 2 * m_Table62.alpha[i] * (delta - m_Table62.epsilon[i])) * (m_Table62.t[i] / tau - 2 * m_Table62.beta[i] * (tau - m_Table62.gamma[i]));
        }
        double psi = 0, theta = 0, Delta = 0, dpsiddelta = 0, dpsidtau=0, dDeltaddelta = 0, dDeltabiddelta = 0,dDeltabidtau=0, dDeltabidtauddelta=0,d2psidtauddelta=0,delta_minus_one_squre = 0, tau_minus_one_squre=0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);
            tau_minus_one_squre = pow(tau - 1, 2.0);
            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*tau_minus_one_squre);

            theta = (1-tau) + m_Table62.A[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i]);

            Delta = pow(theta, 2.0) + m_Table62.B[i]*pow(delta_minus_one_squre, m_Table62.a[i]);
            dpsiddelta = -2*m_Table62.C[i]*(delta - 1) * psi;
            dpsidtau = -2*m_Table62.D[i]*(tau - 1) * psi;
            d2psidtauddelta = 4 * m_Table62.C[i] * psi * m_Table62.D[i] * (delta - 1) * (tau - 1);
            dDeltaddelta = 2*(delta - 1)*(theta * m_Table62.A[i]/m_Table62.beta[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i] - 1) + m_Table62.B[i]*m_Table62.a[i]*pow(delta_minus_one_squre, m_Table62.a[i]-1));
            // if(Delta == 0)
            // {
            //     dDeltabiddelta = 0;
            // }else
            // {
                dDeltabiddelta = dDeltaddelta*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-1);
                dDeltabidtau = -2*theta*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-1);
                dDeltabidtauddelta = -m_Table62.A[i] * m_Table62.b[i] * 2 / m_Table62.beta[i] * pow(Delta, m_Table62.b[i]-1.0) * (delta - 1) * pow(pow(delta - 1, 2.0), 0.5/m_Table62.beta[i] - 1.0) - 2 * theta * m_Table62.b[i] * (m_Table62.b[i] - 1) * pow(Delta, m_Table62.b[i]-2.0) * dDeltaddelta;
            // }
            
            sum4 += m_Table62.n[i]*(pow(Delta,m_Table62.b[i]) * (dpsidtau + delta * d2psidtauddelta) + delta * dDeltabiddelta * dpsidtau + dDeltabidtau * (psi + delta * dpsiddelta) + dDeltabidtauddelta * delta * psi);
        }
        
        return sum1+sum2+sum3+sum4;
    }
    double cH2O::Phi_o(double delta, double tau)
    {
        double sum = 0;
        for (size_t i = 3; i < 8; i++)
        {
            sum += m_Table61.n0[i]*log(1 - exp(-m_Table61.gamma0[i]*tau));
        }
        return log(delta) + m_Table61.n0[0] + m_Table61.n0[1]*tau + m_Table61.n0[2]*log(tau) + sum;
    }
    double cH2O::Phi_o_tau(double delta, double tau)
    {
        double sum = 0;
        for (size_t i = 3; i < 8; i++)
        {
            sum += m_Table61.n0[i]*m_Table61.gamma0[i]*( 1.0/(1.0 - exp(-m_Table61.gamma0[i] * tau)) - 1.0);
        }
        
        return m_Table61.n0[1] + m_Table61.n0[2]/tau + sum;
    }
    double cH2O::Phi_o_tautau(double delta, double tau)
    {
        double sum = 0, exp_gamma_tau=0;
        for (size_t i = 3; i < 8; i++)
        {
            exp_gamma_tau = exp(-m_Table61.gamma0[i] * tau);
            sum += m_Table61.n0[i]*pow(m_Table61.gamma0[i], 2.0)*exp_gamma_tau*pow(1-exp_gamma_tau, -2.0);
        }
        return -m_Table61.n0[2]/pow(tau, 2.0) - sum;
    }
    double cH2O::Phi_r(double delta, double tau)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
        for (size_t i = 0; i < 7; i++)
        {
            sum1 += m_Table62.n[i] * pow(delta, m_Table62.d[i]) * pow(tau, m_Table62.t[i]);
        }
        for (size_t i = 7; i < 51; i++)
        {
            sum2 += m_Table62.n[i] * exp(-pow(delta,m_Table62.c[i])) * 
                    pow(delta, m_Table62.d[i]) * pow(tau, m_Table62.t[i]);
        }
        for (size_t i = 51; i < 54; i++)
        {
            sum3 += m_Table62.n[i] * pow(delta,m_Table62.d[i]) * pow(tau, m_Table62.t[i]) * exp(-m_Table62.alpha[i]*pow(delta - m_Table62.epsilon[i], 2.0) - m_Table62.beta[i]*pow(tau - m_Table62.gamma[i], 2.0) );
        }
        double psi = 0, theta = 0, Delta = 0, dpsidtau = 0, dDeltabidtau = 0, delta_minus_one_squre = 0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);
            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*pow(tau - 1, 2.0));
            theta = (1-tau) + m_Table62.A[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i]);
            Delta = pow(theta, 2.0) + m_Table62.B[i]*pow(delta_minus_one_squre, m_Table62.a[i]);
            sum4 += m_Table62.n[i]*delta*pow(Delta, m_Table62.b[i])*psi;
        }
        
        return sum1+sum2+sum3+sum4;
    }
    double cH2O::Phi_r_tau(double delta, double tau)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
        for (size_t i = 0; i < 7; i++)
        {
            sum1 += m_Table62.n[i] * m_Table62.t[i] * pow(delta, m_Table62.d[i]) * pow(tau, m_Table62.t[i]-1.0);
        }
        for (size_t i = 7; i < 51; i++)
        {
            sum2 += m_Table62.n[i] * m_Table62.t[i] * exp(-pow(delta,m_Table62.c[i])) * 
                    pow(delta, m_Table62.d[i]) * pow(tau, m_Table62.t[i]-1.0);
        }
        for (size_t i = 51; i < 54; i++)
        {
            sum3 += m_Table62.n[i] * pow(delta,m_Table62.d[i]) * pow(tau, m_Table62.t[i]) * exp(-m_Table62.alpha[i]*pow(delta - m_Table62.epsilon[i], 2.0) - m_Table62.beta[i]*pow(tau - m_Table62.gamma[i], 2.0) ) * (m_Table62.t[i]/tau - 2*m_Table62.beta[i]*(tau - m_Table62.gamma[i]) );
        }
        double psi = 0, theta = 0, Delta = 0, dpsidtau = 0, dDeltabidtau = 0, delta_minus_one_squre = 0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);
            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*pow(tau - 1, 2.0));
            theta = (1-tau) + m_Table62.A[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i]);
            Delta = pow(theta, 2.0) + m_Table62.B[i]*pow(delta_minus_one_squre, m_Table62.a[i]);

            dpsidtau = -2*m_Table62.D[i]*(tau - 1) * psi;
            dDeltabidtau = -2*theta*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-1);
            
            sum4 += m_Table62.n[i]*delta*(dDeltabidtau * psi + pow(Delta, m_Table62.b[i])*dpsidtau);
        }
        
        return sum1+sum2+sum3+sum4;
    }
    double cH2O::Phi_r_tautau(double delta, double tau)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
        for (size_t i = 0; i < 7; i++)
        {
            sum1 += m_Table62.n[i] * m_Table62.t[i] * (m_Table62.t[i] - 1.0) * pow(delta, m_Table62.d[i]) * pow(tau, m_Table62.t[i]-2.0);
        }
        for (size_t i = 7; i < 51; i++)
        {
            sum2 += m_Table62.n[i] * m_Table62.t[i] * (m_Table62.t[i] - 1.0) * exp(-pow(delta,m_Table62.c[i])) * 
                    pow(delta, m_Table62.d[i]) * pow(tau, m_Table62.t[i]-2.0);
        }
        for (size_t i = 51; i < 54; i++)
        {
            sum3 += m_Table62.n[i] * pow(delta,m_Table62.d[i]) * pow(tau, m_Table62.t[i]) * exp(-m_Table62.alpha[i]*pow(delta - m_Table62.epsilon[i], 2.0) - m_Table62.beta[i]*pow(tau - m_Table62.gamma[i], 2.0) ) * (pow(m_Table62.t[i]/tau - 2*m_Table62.beta[i]*(tau - m_Table62.gamma[i]), 2.0) - m_Table62.t[i]/pow(tau,2.0) - 2.0*m_Table62.beta[i]);
        }
        double psi = 0, theta = 0, Delta = 0, dpsidtau = 0, dDeltabidtau = 0, d2psidtau2=0, d2Deltabidtau2=0, delta_minus_one_squre = 0, tau_minus_one_squre=0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);
            tau_minus_one_squre = pow(tau - 1, 2.0);
            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*tau_minus_one_squre);
            theta = (1-tau) + m_Table62.A[i]*pow(delta_minus_one_squre, 0.5/m_Table62.beta[i]);
            Delta = pow(theta, 2.0) + m_Table62.B[i]*pow(delta_minus_one_squre, m_Table62.a[i]);

            dpsidtau = -2*m_Table62.D[i]*(tau - 1) * psi;
            dDeltabidtau = -2*theta*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-1);
            d2psidtau2 = 2*m_Table62.D[i]*psi*(2*m_Table62.D[i]*tau_minus_one_squre - 1);
            d2Deltabidtau2 = 2*m_Table62.b[i]*pow(Delta, m_Table62.b[i]-2) *(Delta + 2*theta*theta*(m_Table62.b[i] - 1));
            sum4 += m_Table62.n[i]*delta*(d2Deltabidtau2 * psi + pow(Delta, m_Table62.b[i])*d2psidtau2 + 2*dDeltabidtau*dpsidtau);
        }
        
        return sum1+sum2+sum3+sum4;
    }
    double cH2O::SpecificEnthalpy_T_Rho(double T, double Rho)
    {
        double T_K = T + Kelvin;
        double delta = Rho/Rho_Critic;
        double tau = T_Critic_K/T_K;
        // Note that, R_const in unit of kJ/kg/K
        return R_const*T_K*(1 + tau*(Phi_o_tau(delta, tau) + Phi_r_tau(delta, tau)) + delta*Phi_r_delta(delta, tau) );
    }
    /**
     * \image html water_h.svg "Specific enthalpy of water calculated by swEOS" width=50%.
     * \note The result is compared with python package of <a href="https://iapws.readthedocs.io/en/latest/iapws.iapws95.html#">iapws.IAPWS95</a>
     */
    double cH2O::SpecificEnthalpy(double T, double P)
    {
        return SpecificEnthalpy_T_Rho(T, Rho(T, P));
    }
    double cH2O::Cv_T_Rho(double T, double Rho)
    {
        double T_K = T + Kelvin;
        double delta = Rho/Rho_Critic;
        double tau = T_Critic_K/T_K;
        return (-tau*tau * (Phi_o_tautau(delta, tau) + Phi_r_tautau(delta, tau))) * R_const;
    }
    /**
     * \image html water_cv.svg "Isochoric heat capacity of water calculated by swEOS" width=50%.
     * \note The result is compared with python package of <a href="https://iapws.readthedocs.io/en/latest/iapws.iapws95.html#">iapws.IAPWS95</a>
     */
    double cH2O::Cv(double T, double P)
    {
        return Cv_T_Rho(T, Rho(T, P));
    }
    double cH2O::Cp_T_Rho(double T, double Rho)
    {
        double T_K = T + Kelvin;
        double delta = Rho/Rho_Critic;
        double tau = T_Critic_K/T_K;
        double Cv=Cv_T_Rho(T, Rho);
        return Cv + R_const*pow(1+delta*Phi_r_delta(delta,tau) - delta*tau*Phi_r_deltatau(delta,tau), 2.0)/(1 + 2*delta*Phi_r_delta(delta,tau) + delta*delta*Phi_r_deltadelta(delta,tau));
    }
    /**
     * \image html water_cp.svg "Isobaric heat capacity of water calculated by swEOS" width=50%.
     * \note The result is compared with python package of <a href="https://iapws.readthedocs.io/en/latest/iapws.iapws95.html#">iapws.IAPWS95</a>
     */
    double cH2O::Cp(double T, double P)
    {
        return Cp_T_Rho(T, Rho(T, P));
    }
    double cH2O::mu_T_Rho(double T, double Rho)
    {
        double T_K = T + Kelvin;
        double T_bar=T_K/T_Critic_K;        //T* is critical temperature
        double Rho_bar=Rho/Rho_Critic;      //rho* is the critical density

        // 1.1 mu0
        double mu0=0;
        for (size_t i = 0; i < 4; i++)
        {
            mu0 += m_Table235.H[i]/pow(T_bar,i); //eq. 11
        }
        mu0 = 100*sqrt(T_bar)/mu0; //eq. 11
        // 1.2 mu1
        double mu1=0, mu1_inner=0;
        for (size_t i = 0; i < 6; i++)
        {
            mu1_inner=0;
            for (size_t j = 0; j < 7; j++)
            {
                mu1_inner += m_Table235.Hij[i][j]*pow(Rho_bar - 1, j);
            }
            mu1 += pow(1.0/T_bar - 1.0, i)*mu1_inner;
        }
        mu1 = exp(mu1*Rho_bar);
        // 1.3 mu2
        double Rho2=Rho*0.9995; // small change of rho which is used to calculate gradient(derivative)
        double Tbig_K = T_Critic_K*1.5;//eq. 28
        double Tbig_C = Tbig_K - Kelvin;

        double chi1 = ((Rho-Rho2)/(Pressure_T_Rho(T, Rho)-Pressure_T_Rho(T, Rho2)))*P_Critic/Rho_Critic;
        double chi2 = ((Rho-Rho2)/(Pressure_T_Rho(Tbig_C, Rho)-Pressure_T_Rho(Tbig_C, Rho2)))*P_Critic/Rho_Critic;
        double chi_bar = (chi1 - chi2*Tbig_K/T_K)*Rho_bar; //eq. 28
        if(chi_bar<0)chi_bar=0;

        double xi = m_Table235.xi0 * pow(chi_bar/m_Table235.Gamma0, m_Table235.nu/m_Table235.gamma);
        if(xi<0)xi=0;
        double Y=0;
        const double q_C_xi = m_Table235.q_C*xi;
        const double q_C_xi_square = q_C_xi*q_C_xi;
        const double q_D_xi = m_Table235.q_D*xi;
        const double q_D_xi_square = q_D_xi*q_D_xi;
        if(xi>=0 && xi<=0.3817016416) //page 114 and figure 9
        {
            Y = 0.2 * q_C_xi * pow(q_D_xi, 5.0) * (1 - q_C_xi + q_C_xi_square - 765.0 / 504.0 * q_D_xi_square); //eq 20
        }else
        {
            double psi_D = acos(1.0/sqrt(1 + q_D_xi_square));
            double omega = sqrt(fabs((q_C_xi-1)/(q_C_xi+1)))*tan(psi_D/2.0);
            double L_omega = 0;
            if(q_C_xi>1)
            {
                L_omega = log((1+omega)/(1-omega));
            }else
            {
                L_omega = 2*atan(fabs(omega));
            }
            
            Y = 1.0 / 12.0 * sin(3 * psi_D) - 0.25 / q_C_xi * sin(2 * psi_D) + 1.0 / q_C_xi_square * (1 - 1.25 * q_C_xi_square) * sin(psi_D) - 1.0 / pow(q_C_xi, 3.0) * ((1 - 1.5 * q_C_xi_square) * psi_D - pow(fabs(q_C_xi_square - 1), 1.5) * L_omega); //eq. 26
        }
        double mu2 = exp(m_Table235.chi_mu*Y);

        return mu0*mu1*mu2/1E6;
    }
    /**
     * \image html water_mu.svg "Dynamic viscosity of water calculated by swEOS" width=50%.
     * \note The result is compared with python package of <a href="https://iapws.readthedocs.io/en/latest/iapws.iapws95.html#">iapws.IAPWS95</a>
     */
    double cH2O::mu(double T, double P)
    {
        return mu_T_Rho(T, Rho(T, P));
    }
    /**
     * \image html water_beta.svg "Isotermal compressibility of water calculated by swEOS" width=50%.
     * \note The result is compared with python package of <a href="https://iapws.readthedocs.io/en/latest/iapws.iapws95.html#">iapws.IAPWS95</a>
     */
    double cH2O::beta(double T, double P, double dP)
    {
        double rho = Rho(T, P);
        double P2 = P + dP;
        double P_boil = P_Boiling(T);
        if((sign(P - P_boil)!=sign(P2 - P_boil)) && ((rho - Rho_Critic)<0) )
        {
            P2 = P - dP;
        }
        double rho2 = Rho(T, P2);

        return 1.0/rho * (rho-rho2)/(P-P2) * 1E-5;
    }
    /**
     * \image html water_alpha.svg "Isobaric expansivity of water calculated by swEOS" width=50%.
     * \note The result is compared with python package of <a href="https://iapws.readthedocs.io/en/latest/iapws.iapws95.html#">iapws.IAPWS95</a>
     */
    double cH2O::alpha(double T, double P, double dT)
    {
        double rho = Rho(T, P);
        double T2 = T + dT;
        double P_boil = P_Boiling(T);
        double P2_boil = P_Boiling(T2);
        if ((sign(P - P_boil) != sign(P-P2_boil)) && (P!=P_boil))
        {
            T2 = T - dT;
        }
        double rho2 = Rho(T2, P);
        return -1.0/rho * (rho - rho2)/(T - T2);
    }
} // namespace H2O
