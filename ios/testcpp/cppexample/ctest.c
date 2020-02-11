//
//  ctest.c
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#include "ctest.h"
#include "steam4.h"
void helloC()
{
    printf("hello C\n");
    double T, d, p, s, h, dp, ds, dh;
    Prop *prop0;
    dp = 1.0e-8;
    ds = 1.0e-8;
    dh = 1.0e-8;

    prop0 = newProp('t', 'p', 1);

    T = 273.2;
    d = 0.0;
    p = 1.1 * 1e5;
//    h = 2300.0 * 1e3;
    water_tp(T,p,d,dp,prop0);
//    dumpProp(stdout, prop0);
    printf("rho: %lf, h: %lf mu: %lf\n",prop0->d, prop0->h, viscos(prop0));
    prop0 = freeProp(prop0);
    
    double rho, mu;
    rho_h_mu_Tp(T,p,&rho,&h,&mu);
    printf("func---rho: %lf, h: %lf mu: %lf\n",rho, h, mu);
}

void cal_Tp(double T, double p)
{
    printf("hello C\n");
    double rho, mu, h;
    rho_h_mu_Tp(T,p,&rho,&h,&mu);
    printf("func---rho: %lf, h: %lf mu: %lf\n",rho, h, mu);
}

void rho_h_mu_Tp(double T, double p, double* rho, double* h, double* mu)
{
    Prop *prop0;
    double dp, ds, dh,d;
    dp = 1.0e-8;
    ds = 1.0e-8;
    dh = 1.0e-8;
    d = 0;
    prop0 = newProp('t', 'p', 1);
    
    water_tp(T,p,d,dp,prop0);
    *rho=prop0->d;
    *h=prop0->h;
    *mu=viscos(prop0);
}

