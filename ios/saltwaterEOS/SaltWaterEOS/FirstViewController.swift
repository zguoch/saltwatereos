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
        
        let quote = "Haters gonna hate"
        let font = UIFont.systemFont(ofSize: 12)
        let attributes: [NSAttributedString.Key: Any] = [
            .font: font,
            .foregroundColor: UIColor.red,
        ]
        let attributedQuote = NSAttributedString(string: quote, attributes: attributes)
            
        var value = 50.0
        var x = [1.0,2.0,3.0]
        struct PROPSWEOS{
           var Region: CInt
           var Region_str: NSString
           var Rho: Double
            var Rho_l: Double
            var Rho_v: Double
            var Rho_h: Double
        }
        var prop = PROPSWEOS(Region:0, Region_str: "aa", Rho:0.1, Rho_l: 0, Rho_v: 0, Rho_h: 0)
//        swEOS_wrapper().hellocpp_wrappe(4.4, 5,x);
        swEOS_wrapper().prop_pTX(p, T, 0.032,&prop.Region,&prop.Rho,&prop.Rho_l ,&prop.Rho_v,&prop.Rho_h)
        print("value: ", prop.Rho)
        textView_Results.attributedText=attributedQuote
        textView_Results.text=textView_Results.text+String(format:"\nDensity: %.1f kg/m3\nViscosity: %.1f Pa s\n",prop.Rho,prop.Rho_v);
    }
    
    @IBOutlet weak var textView_Results: UITextView!
    @IBAction func SelectCalculationMode_Scatter(_ sender: Any) {
    }
    
}

