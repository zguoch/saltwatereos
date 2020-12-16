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
        double T_square = T*T;
        double P_square = P*P;
        double T_minus_Ttriple = T - T_Triple;
        double c[3]={226713, 44.6652, -7.41999E-05};
        double r[5]={1148.81, 0.275774, 8.8103E-05,
                    -0.0017099 - 3.82734E-06 * T / 2.0 - 8.65455E-09 * T_square / 3.0,
                    5.29063E-08 - 9.63084E-11 * T / 2.0 + 6.50745E-13 * T_square / 3.0
                    };
        return c[0] + r[0]*T_minus_Ttriple 
                    + r[1]*pow(T_minus_Ttriple, 2.0)
                    + r[2]*pow(T_minus_Ttriple, 3.0)
                    + (c[1] + r[3]*T)*P 
                    + (c[2] + r[4]*T)*P_square;
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
