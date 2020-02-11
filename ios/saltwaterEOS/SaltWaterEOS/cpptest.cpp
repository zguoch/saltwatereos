//
//  cpptest.cpp
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#include "cpptest.hpp"
//#include "steam4.h"
#include "ctest.h"
#include "H2ONaCl.h"

void hellocpp()
{
    cout<<"hello cpp\n";
    
    cout.precision(8);//control cout precision of float
    cH2ONaCl eos(31600000,100,0.3);
    eos.Calculate();
    
    cout<<"Region: "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
    cout<<"Rho_l: "<<eos.m_prop.Rho_l<<" Rho_v: "<<eos.m_prop.Rho_v<<" Rho_h: "<<eos.m_prop.Rho_h<<endl;
    cout<<"h_l: "<<eos.m_prop.H_l<<" h_v: "<<eos.m_prop.H_v<<" h_h: "<<eos.m_prop.H_h<<endl;
    cout<<"Rho: "<<eos.m_prop.Rho<<" H: "<<eos.m_prop.H<<endl;
    cout<<"S_l: "<<eos.m_prop.S_l<<" S_v: "<<eos.m_prop.S_v<<" S_h: "<<eos.m_prop.S_h<<endl;
    cout<<"Mu_l: "<<eos.m_prop.Mu_l<<" Mu_v: "<<eos.m_prop.Mu_v<<endl;
    
    helloC();
//    double T, d, p, s, h, dp, ds, dh;
//    Prop *prop0;
//    dp = 1.0e-8;
//    ds = 1.0e-8;
//    dh = 1.0e-8;
//
//    prop0 = newProp('t', 'p', 1);
    
//    T = 273.2;
//    d = 0.0;
//    p = 1.1 * 1e5;
//    //    h = 2300.0 * 1e3;
//    water_tp(T,p,d,dp,prop0);
//    //    dumpProp(stdout, prop0);
//    printf("rho: %lf, h: %lf mu: %lf\n",prop0->d, prop0->h, viscos(prop0));
//    prop0 = freeProp(prop0);
    
}
