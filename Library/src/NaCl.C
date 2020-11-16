#include "NaCl.H"

namespace NaCl
{
    cNaCl::cNaCl(/* args */)
    {
    }
    
    cNaCl::~cNaCl()
    {
    }
    double cNaCl::P_Boiling(double T)
    {
        return pow(10, (log10(P_Triple) + 9418.12 * (1.0 / (T_Triple + 273.15) - 1.0 / (T + 273.15))));
    }
} // namespace NaCl
