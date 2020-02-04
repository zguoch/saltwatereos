//
//  DataViewController.m
//  callCpp
//
//  Created by Zhikui on 2/3/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

#import "DataViewController.h"
#import "testcpp.hpp"

@interface DataViewController ()
{
    Greeting greeting;
}

@end

@implementation DataViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view.
}


- (void)viewWillAppear:(BOOL)animated {
    [super viewWillAppear:animated];
    self.dataLabel.text = [self.dataObject description];
}


@end
