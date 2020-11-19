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
     * 
     */
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
        double psi = 0, theta = 0, Delta = 0, dpsiddelta = 0, dDeltaddelta = 0, dDeltabiddelta = 0, delta_minus_one_squre = 0;
        for (size_t i = 54; i < 56; i++)
        {
            delta_minus_one_squre = pow(delta - 1, 2.0);

            psi = exp(-m_Table62.C[i]*delta_minus_one_squre - m_Table62.D[i]*pow(tau - 1, 2.0));

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
    double cH2O::Rho(double T, double P)
    {//DEBUG: Rho_Water# in water_prop.vb
        double T_K = T + Kelvin;
        double Rho1 = 0, Rho2 = 0;
        if(T_K <= T_Critic_K)
        {
            double RoundDown_P = floor(P*1000)/1000.0;//DEBUG
            double RoundDown_BoilingCurve = floor(BoilingCurve(T)*1000)/1000.0;//DEBUG
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
} // namespace H2O
