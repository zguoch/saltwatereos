#include "H2ONaCl.H"
#include <iostream>
int main()
{
    H2ONaCl::cH2ONaCl eos;
    double p=200E5;   //Pa
    double T=400+273.15;   //K
    double X=0.2;    //wt.% NaCl
    double rho=eos.rho_pTX(p,T,X); //kg/m3
    std::cout<<" Pressure(bar): "<<p/1E5<<"\n"
             <<" Temperature(deg.C): "<<T-273.15<<"\n"
             <<" Salinity (wt.% NaCl): "<<X<<std::endl;
    std::cout<<" Density(kg/m3): "<<rho<<std::endl;
}