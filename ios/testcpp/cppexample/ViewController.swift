//
//  ViewController.swift
//  cppexample
//
//  Created by Zhikui on 2/4/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

import UIKit

class ViewController: UIViewController {

    override func viewDidLoad() {
        super.viewDidLoad()
        // Do any additional setup after loading the view.
        
//        helloC();
        cpp_wrapper().hellocpp_wrappe();
        
//        let T:Double=375;
//        let p:Double=100e5;
//        var rho:Double=0;
//        var h:Double=0;
//        var mu:Double=0;
////        cal_Tp(T, p);
//        rho_h_mu_Tp(T,p,&rho,&h,&mu);
//        print("T: ", T, "p: ",p, "rho: ",rho, "h: ", h, "mu: ", mu);
    }
    @IBOutlet weak var Tbox_rho: UITextField!
    @IBAction func calculate(_ sender: Any) {
        
                let T:Double=375;
                let p:Double=100e5;
                var rho:Double=0;
                var h:Double=0;
                var mu:Double=0;
        //        cal_Tp(T, p);
                rho_h_mu_Tp(T,p,&rho,&h,&mu);
                print("T: ", T, "p: ",p, "rho: ",rho, "h: ", h, "mu: ", mu);
        
        Tbox_rho.text=String(format:"%.1f",rho);
//        let alert = UIAlertController(title: "Did you bring your towel?", message: "It's recommended you bring your towel before continuing.", preferredStyle: .alert)
//
//        alert.addAction(UIAlertAction(title: "Yes", style: .default, handler: nil))
//        alert.addAction(UIAlertAction(title: "No", style: .cancel, handler: nil))
//
//        self.present(alert, animated: true)
    }
    
}

