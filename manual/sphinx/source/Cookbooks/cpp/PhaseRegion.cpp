#include "H2ONaCl.H"
#include <iostream>
H2ONaCl::cH2ONaCl eos;

int main()
{
    H2ONaCl::cH2ONaCl eos;
    double p=200; //bar
    double T=400; //deg.C
    double X=0.032; //wt.% NaCl
    H2ONaCl::PhaseRegion region=eos.findPhaseRegion(T, p, X);
    std::string region_name=eos.getPhaseRegionName(region);
    std::cout<<" Pressure(bar): "<<p<<"\n"
             <<" Temperature(deg.C): "<<T<<"\n"
             <<" Salinity (wt.% NaCl): "<<X<<std::endl;
    std::cout<<" Phase region index: "<<region<<std::endl;
    std::cout<<" Phase region name: "<<region_name<<std::endl;
}