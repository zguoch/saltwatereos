//
//  ctest.h
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#ifndef ctest_h
#define ctest_h

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif
    
    extern void helloC(void);
    
#ifdef __cplusplus
} // extern "C"
#endif


void rho_h_mu_Tp(double T, double p, double* rho, double* h, double* mu);
void cal_Tp(double T, double p);
#endif /* ctest_h */
