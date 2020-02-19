//
//  cpp_wrapper.h
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface swEOS_wrapper : NSObject

-(double ) hellocpp_wrappe:(double)a:(double)b:(double*)Rho;

-(void) prop_pTX:(double)P:(double)T:(double)X:(double*)Rho;

@end
