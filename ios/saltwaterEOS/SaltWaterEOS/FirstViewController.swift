//
//  FirstViewController.swift
//  SaltWaterEOS
//
//  Created by Zhikui on 2/11/20.
//  Copyright © 2020 Zhikui. All rights reserved.
//

import UIKit

class FirstViewController: UIViewController, UIPickerViewDataSource,UIPickerViewDelegate {
    
    var independentVar=["Temperature (C)","Pressure (bar)", "Salinity"]
    var picker_independentVar=UIPickerView()
    
    
    func numberOfComponents(in pickerView: UIPickerView) -> Int {
        return 1
    }
    
    func pickerView(_ pickerView: UIPickerView, numberOfRowsInComponent component: Int) -> Int {
        return independentVar.count
    }
    func pickerView(_ pickerView: UIPickerView, didSelectRow row: Int, inComponent component: Int) {
//        textField_independentVars.text=independentVar[row]
    }
    func pickerView(_ pickerView: UIPickerView, titleForRow row: Int, forComponent component: Int) -> String? {
        return independentVar[row]
    }

    override func viewDidLoad() {
        super.viewDidLoad()
        // Do any additional setup after loading the view.
        picker_independentVar.delegate=self
        picker_independentVar.dataSource=self
    }

    @IBOutlet weak var textField_ThirdVar: UITextField!
    @IBOutlet weak var textField_SecondVar: UITextField!
    @IBOutlet weak var textField_FirstVar: UITextField!
    @IBAction func Calculate(_ sender: Any) {
        let T = NSString(string: textField_SecondVar.text!).doubleValue;
        let p = NSString(string: textField_FirstVar.text!).doubleValue;
                var rho:Double=0;
                var h:Double=0;
                var mu:Double=0;
        //        cal_Tp(T, p);
                rho_h_mu_Tp(T,p*1e5,&rho,&h,&mu);
                print("T: ", T, "p: ",p, "rho: ",rho, "h: ", h, "mu: ", mu);
        textView_Results.text=textView_Results.text+String(format:"\nDensity: %.1f kg/m3\nViscosity: %.1f Pa s\n",rho,rho);
    }
    
    @IBOutlet weak var textView_Results: UITextView!
    @IBAction func SelectCalculationMode_Scatter(_ sender: Any) {
    }
    
}

