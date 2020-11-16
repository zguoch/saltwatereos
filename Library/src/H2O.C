#include "H2O.H"

namespace H2O
{
    cH2O::cH2O(/* args */)
    {
    }
    
    cH2O::~cH2O()
    {
    }
    double cH2O::P_Boiling(double T)
    {
        double T_K =T;
        if(T_K==0)T_K=0.01;
        T_K = T_K + Kelvin;
        double T_inv = 1 - T_K / 647.096;
        double a[6] = {-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502};
        double P_boil = 0;

        P_boil = exp((T_Critic + Kelvin) / T_K * (a[0] * T_inv + a[1] * pow(T_inv,1.5) + a[2] * pow(T_inv, 3.0) + a[3] * pow(T_inv, 3.5) + a[4] * pow(T_inv,4.0) + a[5] * pow(T_inv ,7.5))) * P_Critic;
        return P_boil;
    }
} // namespace H2O
