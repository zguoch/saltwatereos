#include "NaCl.H"

namespace NaCl
{
    cNaCl::cNaCl(/* args */)
    {
    }
    
    cNaCl::~cNaCl()
    {
    }
    /**
     * \f$ log_{10}(P_{NaCl, liquid}) = log_{10}(P_{triple, liquid}) + b_{boil}\left( \frac{1}{T_{triple, NaCl} + 273.15} - \frac{1}{T + 273.15} \right) , b_{boil} = 0.941812\times 10^4 \f$
     * 
     * \image html Driesner_Heinrich_Fig4.png "Vapor pressure curves (sublimation and boiling curves) of NaCl." width=25%.
     * Vapor pressure curves (sublimation and boiling curves) of NaCl, figure 4 of reference \cite Driesner2007Part1. 
     * 
     * \image html HaliteSublimationBoilingCurves.svg "Vapor pressure curves (sublimation and boiling curves) of NaCl." width=25%.
     * see also #P_Sublimation 
     */
    double cNaCl::P_Boiling(double T)
    {
        return pow(10, (log10(P_Triple) + 9418.12 * (1.0 / (T_Triple + 273.15) - 1.0 / (T + 273.15))));
    }
    /**
     * \f{equation}
     * \rho_{halite} = \rho_{halite}^0 + l(T)P
     * \f}
     * where \f$ \rho_{halite}^0 \f$ and \f$ l(T) \f$, and corresponding coefficient can be found in equation (2), (3) and Table 3 of \cite Driesner2007Part2. 
     */
    double cNaCl::Rho_Solid(double T, double P)
    {
        double l[6] = {2170.4, -0.24599, -9.5797E-05, 0.005727, 0.002715, 733.4};
        double rho0 = l[0] + l[1]*T +l[2]*T*T;
        double ll = l[3] + l[4]*exp(T/l[5]);
        return rho0 + ll*P;
    }
    /**
     * \f{equation}
     * \rho_{NaCl, liquid} = \frac{\rho_{NaCl, liquid}^0}{1 - 0.1ln(1+10P \kappa_{NaCl, liquid})}
     * \f}
     * where \f$ \rho_{NaCl, liquid}^0 \f$, \f$\kappa_{NaCl, liquid} \f$ (compressibility) and \f$ m_i(i=0,...,5) \f$ can be found in equation (5), (6) and Table 3 of \cite Driesner2007Part2, respectively. 
     */
    double cNaCl::Rho_Liquid(double T, double P)
    {
        double m[6] = {58443, 23.772, 0.018639, -1.9687E-06, -1.5259E-05, 5.5058E-08};
        double rho0 = m[0]/(m[1] + m[2]*T + m[3]*T*T);
        double kappa = m[4] + m[5]*T;
        return rho0/(1 - 0.1*log(1 + 10*P*kappa));
    }
    double cNaCl::SpecificEnthalpy(double T, double P)
    {
        double l0  = 2.1704e3;
        double l1  = -2.4599e-1;
        double l2  = -9.5797e-5;
        double l3  = 5.727e-3;
        double l4  = 2.715e-3;
        double l5  = 733.4;
        //JH 02/2021: bugfix (must be T^3 in calculating F_Cp, not T^2)
        double T_100      = 100;
        double rho0_H     = l0 + l1*T_100 + l2*T_100*T_100;
        double drho0_dT_H = l1 + 2 * l2 * T_100; // % d/dT(rho_H) at T=100 C
        double l          = l3 + l4 * exp(T_100/l5); // Eq.(3)
        double dl_dT      = l4/l5 * exp(T_100/l5); // dl/dT at T=100 C
        // P_100      = 100;
        // F_V_100    = 1./(l * 1e-5) .* log(rho0_H  +  l .* P_100); // 1e-5 comes from integration for P in Pa and not in bar
        double F_V_100    = 86792219.2899765; // result of the above commented lines
        // P_100      = 100;
        // dF_V_100dT =  1e5 * ( ( - dl_dT ./ (l^2) ) .* log(rho0_H  +  l .* P_100) ...
        //               +  1./l .* ( drho0_dT_H + dl_dT .* P_100 ) ./ (rho0_H  +  l .* P_100) ); 
        double dF_V_100dT = -43058.0130278123; // result of the above commented lines
        double F_V_p      = 1/(l* 1e-5) * log(rho0_H  +  l * P); 
        double dF_V_pdT   = 1e5 * ( ( - dl_dT / (l*l) ) * log(rho0_H  +  l * P) + 1/(l) * ( drho0_dT_H + dl_dT * P ) / (rho0_H  +  l * P) );       
        double h_H_ref    = 9.415867359e4; // from Driesner's table (which one ???)
        double h_H_p_t100 = h_H_ref + (F_V_p - F_V_100 - (T_100 + Kelvin) * (dF_V_pdT - dF_V_100dT));
        // Table 5
        double r0      = 1148.81;
        double r1      = 0.275774;
        double r2      = 8.8103e-5;
        double r3a     = -1.7099e-3;
        double r3b     = -3.82734e-6;
        double r3c     = -8.65455e-9;
        // r3      = r3a  + r3b .* T  + r3c .* T.^2;
        double r4      = 5.29063e-8 - 9.63084e-11 + 6.50745e-13;
        // Equation (30):
        //     Cp_h = r0 + 2.*r1.*(T - T_triple) + 3.*r2.*(T- T_triple).^2 ...
        //             + r3.*P + r4.*(P.^2);

        // ...part 1 (independent of pressure P)
        // F_Cp_100_const = r0 * T_100 + r1.* T_100.^2 - 2 * r1.* T_100 .* T_triple ...
        //                + r2.*T_100.^3 - 3 *r2 .* T_100.^2 .* T_triple + 3 * r2 .* T_100 .* T_triple.^2;
        double F_Cp_100_const = 88393.4640361410; // verified (JH 02/2021)

        //...part 2 (depending on P)
        // value of antiderivative of (30) at (T==100 C,P)

        // r3_term = r3a .* T_100 + 0.5 * r3b .* T_100^2 + 0.5 * r3c.*T_100^2; WRONG
        // r3_term = -0.190169972750000; // WRONG VALUE

        // r3_term = r3a .* T_100 + 0.5 * r3b .* T_100^2 + (1/3) * r3c .* T_100^3; 
        double r3_term = -0.193011550000000; // verified (JH 02/2021)

        // r4_term = r4 .* T_100;
        double r4_term = 5.28106423450000e-6; // verified (JH 02/2021)

        // putting it together
        double F_Cp_T100 = F_Cp_100_const + r3_term * P  + r4_term * P*P;
        double T_square = T*T;
        double F_Cp_at_T = r0 * T + r1 * T_square - 2 * r1 * T  * T_Triple 
          + r2 *T_square*T - 3 *r2  * T_square  * T_Triple + 3 * r2  * T  * T_Triple*T_Triple
          + (r3a  * T + 0.5 * r3b  * T_square + (1/3.0) * r3c  * T_square*T)  * P 
          + r4  * T  * P*P;

        return h_H_p_t100 + F_Cp_at_T - F_Cp_T100;
    }
    double cNaCl::Cp(double T, double P)
    {
        double T_square = T*T;
        double P_square = P*P;
        double T_minus_Ttriple = T - T_Triple;
        double r[5] = {1148.81, 0.275774, 8.8103E-05, 
                    -0.0017099 - 3.82734E-06 * T - 8.65455E-09 * T_square,
                    5.29063E-08 - 9.63084E-11 * T + 6.50745E-13 * T_square};
        return r[0] 
               + 2 * r[1] * T_minus_Ttriple 
               + 3 * r[2] * T_minus_Ttriple*T_minus_Ttriple 
               + r[3] * P 
               + r[4] * P_square;
    }
} // namespace NaCl
