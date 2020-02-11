//
//  FirstViewController.swift
//  SaltWaterEOS
//
//  Created by Zhikui on 2/11/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

import UIKit

class FirstViewController: UIViewController {

    override func viewDidLoad() {
        super.viewDidLoad()
        // Do any additional setup after loading the view.
    }

    @IBOutlet weak var testRhoText: UITextField!
    @IBAction func Calculate(_ sender: Any) {
        let T:Double=375;
                let p:Double=100e5;
                var rho:Double=0;
                var h:Double=0;
                var mu:Double=0;
        //        cal_Tp(T, p);
                rho_h_mu_Tp(T,p,&rho,&h,&mu);
                print("T: ", T, "p: ",p, "rho: ",rho, "h: ", h, "mu: ", mu);
        
        testRhoText.text=String(format:"%.1f",rho);
    }
    
}

