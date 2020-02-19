//
//  cpp_wrapper.h
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#import <Foundation/Foundation.h>


@interface swEOS_wrapper : NSObject

-(void) prop_pTX:(double)P:(double)T:(double)X:(int*)Region:(double*)Rho:(double*)Rho_l:(double*)Rho_v:(double*)Rho_h;

@end
