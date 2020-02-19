//
//  cpp_wrapper.m
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#import "swEOS_wrapper.h"
#import "H2ONaCl.h"

@implementation swEOS_wrapper

-(void) prop_pTX:(double)P:(double)T:(double)X:(int*)Region:(double*)Rho:(double*)Rho_l:(double*)Rho_v:(double*)Rho_h
{
    SWEOS::cH2ONaCl eos;
    eos.prop_pTX(P*1e5, T, X);
    *Rho=eos.m_prop.Rho;
    *Rho_v=eos.m_prop.Rho_v;
    cout<<"test print"<<endl;
    cout<<"P: "<<P<<" T: "<<T<<" X: "<<X<<" Rho: "<<eos.m_prop.H<<" Region: "<<eos.m_prop.Region<<endl;
    
}

@end
