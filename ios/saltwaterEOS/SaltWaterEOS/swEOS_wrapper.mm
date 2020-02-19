//
//  cpp_wrapper.m
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#import "swEOS_wrapper.h"
#import "cpptest.hpp"
#import "H2ONaCl.h"

@implementation swEOS_wrapper

-(double) hellocpp_wrappe:(double)a:(double)b:(double*)Rho
{
    PROP_H2ONaCl prop0=hellocpp(a);
    *Rho=prop0.Rho;
    return 33;
}

-(void) prop_pTX:(double)P:(double)T:(double)X:(double*)Rho
{
    cH2ONaCl eos(31600000,100,0.3);
    eos.Calculate();
}

@end
