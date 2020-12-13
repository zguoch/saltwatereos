#include "H2ONaCl.H"
#include "omp.h"
H2ONaCl::cH2ONaCl eos;
// 1. Halite melting curve, Fig. 3 of Driesner and Heinrich(2007)
void test_HaliteMelting();
// 2. Vapor pressure curves (sublimation and boling curves), Fig. 4 of Driesner and Heinrich(2007)
void test_SublimationBoiling();
// 3. Critical pressure and salinity , Fig. 5 of Driesner and Heinrich(2007)
void test_CriticalPressure_Composition();
// 4. Halite liquidus, Fig. 7 of Driesner and Heinrich(2007)
void test_HaliteLiquidus();
// 5. Halite-saturated vapor composition, Fig. 8 of Driesner and Heinrich(2007)
void test_HaliteSaturatedVaporComposition();
// 6. Pressure at Vapor + Liquid + Halite coexistence. Fig. 9 of Driesner and Heinrich(2007)
void test_P_VLH();
// 7. Salinity on liquid branch of vapor + liquid coexist surface. Fig. 12 of Driesner and Heinrich(2007)
void test_Salinity_VaporLiquidCoexist_LiquidBranch();
// Water cH2O
void test_water_Curves();
void test_water_props(double Tmin, double Tmax, double Pmin, double Pmax, double dP=1E5, double dT=1, bool writeTP = false);
// 8. Volume fraction, Fig. 2 of Driesner(2007)
void testT_V_star();
// test V_extrapol for low TP region, Eq. 17 of Driesner (2007)
void test_V_brine_NaCl_lowThighT();
void test_V_extrapol();
void test_NaClH2O_props(double Tmin, double Tmax, double Pmin, double Pmax, double dP=1E5, double dT=1, bool writeTP = false);

void test_NaClH2O_props(double Tmin, double Tmax, double Pmin, double Pmax, double dP, double dT, bool writeTP)
{
  H2ONaCl::cH2ONaCl eos;
  vector<double> arrX={10, 20};
  for (size_t i = 0; i < arrX.size(); i++)
  {
    double X = arrX[i]/100.0;
    string filename="H2ONaCl_rho_X"+to_string((int)arrX[i])+".dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double P = Pmin; P <= Pmax; P=P+dP)
    {
      for (double T = Tmin; T <= Tmax; T=T+dT)
      {
        double rho = eos.Rho(T, P, X);
        fout<<rho<<" ";
      }
      fout<<"\n";
    }
    fout.close();
  }
  
  if(writeTP)
  {
    string filename_TT="H2ONaCl_T.dat";
    ofstream fout_TT(filename_TT);
    if(!fout_TT)cout<<"Open file failed: "<<filename_TT<<endl;
    string filename_PP="H2ONaCl_P.dat";
    ofstream fout_PP(filename_PP);
    if(!fout_PP)cout<<"Open file failed: "<<filename_PP<<endl;
    for (double T = Tmin; T <= Tmax; T=T+dT)
    {
      fout_TT<<T<<"\n";
    }
    fout_TT.close();
    for (double P = Pmin; P <= Pmax; P=P+dP)
    {
        fout_PP<<P<<"\n";
    }
    fout_PP.close();
  }
}
void test_water_props(double Tmin, double Tmax, double Pmin, double Pmax, double dP, double dT, bool writeTP)
{
  H2ONaCl::cH2ONaCl eos;
  // full range
  {
    string filename="water_rho.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double P = Pmin; P <= Pmax; P=P+dP)
    {
      for (double T = Tmin; T <= Tmax; T=T+dT)
      {
        double rho = eos.m_water.Rho(T, P);
        fout<<rho<<" ";
      }
      fout<<"\n";
    }
    fout.close();
    if(writeTP)
    {
      string filename_TT="water_T.dat";
      ofstream fout_TT(filename_TT);
      if(!fout_TT)cout<<"Open file failed: "<<filename_TT<<endl;
      string filename_PP="water_P.dat";
      ofstream fout_PP(filename_PP);
      if(!fout_PP)cout<<"Open file failed: "<<filename_PP<<endl;
      for (double T = Tmin; T <= Tmax; T=T+dT)
      {
        fout_TT<<T<<"\n";
      }
      fout_TT.close();
      for (double P = Pmin; P <= Pmax; P=P+dP)
      {
          fout_PP<<P<<"\n";
      }
      fout_PP.close();
    }
  }
}
void test_water_Curves()
{
  H2ONaCl::cH2ONaCl eos;
  // boiling curve
  {
    string filename="Water_boiling.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = H2O::T_Triple; T <= H2O::T_Critic; T=T+1)
    {
      // double P=eos.m_water.P_Boiling(T);
      double P=eos.m_water.BoilingCurve(T);
      // double rho_v_sat = eos.m_water.Rho_Vapor_Saturated(T);
      // double rho_l_sat = eos.m_water.Rho_Liquid_Saturated(T);
      // double rho = eos.m_water.Rho(T,P2);
      // double P_T_rho = eos.m_water.Pressure_T_Rho(T, rho_v_sat);
      // fout<<T<<" "<<P<<" "<<P2<<" "<<rho_v_sat<<" "<<rho_l_sat<<" "<<rho<<" "<<P_T_rho<<endl;
      fout<<T<<" "<<P<<endl;
    }
    fout<<H2O::T_Critic<<" "<<H2O::P_Critic<<endl;
    fout.close();
  }
  // sublimation curve
  {
    string filename="Water_sublimation.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = -50; T <= H2O::T_Triple; T=T+0.001)
    {
      double P=eos.m_water.SublimationCurve(T);
      fout<<T<<" "<<P<<endl;
    }
    fout.close();
  }
  // melting curve: ice I
  {
    string filename="Water_iceI.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    // iceI
    for (double T = H2O::T_K_ice_min[H2O::iceI]; T <= H2O::T_Triple_K; T=T+0.1)
    {
      double P=eos.m_water.MeltingCurve(T-Kelvin, true);
      fout<<T-Kelvin<<" "<<P<<endl;
    }
    fout<<H2O::T_Triple<<" "<<H2O::P_Triple<<endl;
    fout.close();
  }
  // melting curve: ice III
  {
    string filename="Water_iceIII.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    double Tend = H2O::T_K_ice_max[H2O::iceIII];
    for (double T = H2O::T_K_ice_min[H2O::iceIII]; T < Tend; T=T+0.01)
    {
      double P=eos.m_water.MeltingCurve(T-Kelvin);
      fout<<T-Kelvin<<" "<<P<<endl;
    }
    fout<<Tend-Kelvin<<" "<<eos.m_water.MeltingCurve(Tend-Kelvin)<<endl;
    fout.close();
  }
  // melting curve: ice III
  {
    string filename="Water_iceV.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    double Tend = H2O::T_K_ice_max[H2O::iceV];
    for (double T = H2O::T_K_ice_min[H2O::iceV]; T < Tend; T=T+0.01)
    {
      double P=eos.m_water.MeltingCurve(T-Kelvin);
      fout<<T-Kelvin<<" "<<P<<endl;
    }
    // fout<<Tend<<" "<<eos.m_water.MeltingCurve(Tend-Kelvin)<<endl;
    fout.close();
  }
  // melting curve: ice VI
  {
    string filename="Water_iceVI.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    double Tend = H2O::T_K_ice_max[H2O::iceVI];
    for (double T = H2O::T_K_ice_min[H2O::iceVI]; T < Tend; T=T+0.1)
    {
      double P=eos.m_water.MeltingCurve(T-Kelvin);
      fout<<T-Kelvin<<" "<<P<<endl;
    }
    // fout<<Tend<<" "<<eos.m_water.MeltingCurve(Tend-Kelvin)<<endl;
    fout.close();
  }
  // melting curve: ice VI
  {
    string filename="Water_iceVII.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    double Tend = H2O::T_K_ice_max[H2O::iceVII];
    for (double T = H2O::T_K_ice_min[H2O::iceVII]; T < Tend; T=T+0.1)
    {
      double P=eos.m_water.MeltingCurve(T-Kelvin);
      fout<<T-Kelvin<<" "<<P<<endl;
    }
    // fout<<Tend<<" "<<eos.m_water.MeltingCurve(Tend-Kelvin)<<endl;
    fout.close();
  }
}
void test_P_VLH()
{
  H2ONaCl::cH2ONaCl eos;
  // full range
  {
    string filename="P_VLH_fullrange.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = 0; T <= 800; T=T+5)
    {
      double P_VLH;
      P_VLH=eos.P_VaporLiquidHaliteCoexist(T);
      fout<<T<<" "<<P_VLH<<endl;
    }
    fout.close();
  }
  
  // close to critical temperature of water
  {
    string filename="P_VLH_lowT.dat";
    ofstream fout(filename); 
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = 0; T <= 300; T=T+0.5)
    {
      double P_VLH;
      P_VLH=eos.P_VaporLiquidHaliteCoexist(T);
      fout<<T<<" "<<P_VLH<<endl;
    }
    fout.close();
  }
}

void test_HaliteMelting()
{
  H2ONaCl::cH2ONaCl eos;
  string filename="HaliteMeltingCurve.dat";
  ofstream fout(filename);
  if(!fout)
  {
    cout<<"Open file failed: "<<filename<<endl;
  }
  for (double P = 0; P <= 5000; P=P+50)
  {
    double T_hm = eos.m_NaCl.T_Melting(P);
    fout<<T_hm<<" "<<P<<endl;
  }
  fout.close();

}
void test_CriticalPressure_Composition()
{
  H2ONaCl::cH2ONaCl eos;
  // full range
  {
    string filename="HaliteCritical_P_X_fullrange.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = H2O::T_Critic; T <= 1000; T=T+10)
    {
      double P,X;
      eos.P_X_Critical(T,P,X);
      fout<<T<<" "<<P<<" "<<X<<endl;
    }
    fout.close();
  }
  
  // close to critical temperature of water
  {
    string filename="HaliteCritical_P_X_closeWaterCriticalT.dat";
    ofstream fout(filename); 
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = H2O::T_Critic; T <= 460; T=T+10)
    {
      double P,X;
      eos.P_X_Critical(T,P,X);
      fout<<T<<" "<<P<<" "<<X<<endl;
    }
    fout.close();
  }
}

void test_SublimationBoiling()
{
  H2ONaCl::cH2ONaCl eos;
  string filename="HaliteSublimationBoilingCurve.dat";
  ofstream fout(filename);
  if(!fout)cout<<"Open file failed: "<<filename<<endl;
  
  for (double T = 300; T <= 1100; T=T+10)
  {
    double P_subl = eos.m_NaCl.P_Sublimation(T);
    double P_boil = eos.m_NaCl.P_Boiling(T);
    fout<<T<<" "<<P_subl<<" "<<P_boil<<endl;
  }
  fout.close();
}

void test_HaliteLiquidus()
{
  std::cout<<"Test halite liquids start"<<endl;
  H2ONaCl::cH2ONaCl eos;
  vector<double> arryP={5, 500, 2000, 4000}; //bar
  for (size_t i = 0; i < arryP.size(); i++)
  {
    string filename="X_HaliteLiquidus_P"+to_string((int)arryP[i])+"bar.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = 5; T < 900; T=T+5) //deg.C
    {
      double X_liquidus=eos.X_HaliteLiquidus(T,arryP[i]); //mol fraction
      fout<<T<<" "<<X_liquidus<<endl;
    }
    fout.close();
  }
  vector<double> arryT={25};
  for (size_t i = 0; i < arryT.size(); i++)
  {
    string filename="X_HaliteLiquidus_T"+to_string((int)arryT[i])+"C.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double P = 0; P < 5000; P=P+100) //deg.C
    {
      double X_liquidus=eos.X_HaliteLiquidus(arryT[i],P); //mol fraction
      fout<<P<<" "<<X_liquidus<<endl;
    }
    fout.close();
  }
  
  
  std::cout<<"Test halite liquids end\n"<<endl;
}

void test_HaliteSaturatedVaporComposition()
{
  std::cout<<"Test halite satureated vapor composition start"<<endl;
  H2ONaCl::cH2ONaCl eos;
  vector<double> arryP={1, 10, 50, 100, 200,300}; //bar
  for (size_t i = 0; i < arryP.size(); i++)
  {
    string filename="X_HaliteSaturatedVapor_P"+to_string((int)arryP[i])+"bar.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = 100; T <= 820; T=T+0.5) //deg.C
    {
      double X_liquidus=eos.X_VaporHaliteCoexist(T,arryP[i]); //mol fraction
      // calculate Halite + Liquid + Vapor coexist composition
      double P_VLH = eos.P_VaporLiquidHaliteCoexist(T);
      double X_VLH = eos.X_VaporHaliteCoexist(T,P_VLH);
      if(X_liquidus>X_VLH)continue;
      fout<<T<<" "<<X_liquidus<<endl;
    }
    fout.close();
  }
  // 
  {
    string filename="X_VLH.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = 100; T <= 800; T=T+0.5) //deg.C
    {
      // calculate Halite + Liquid + Vapor coexist composition
      double P_VLH = eos.P_VaporLiquidHaliteCoexist(T);
      double X_VLH = eos.X_VaporHaliteCoexist(T,P_VLH);
      fout<<T<<" "<<X_VLH<<endl;
    }
    for (double T = 800; T <= 820; T=T+0.01) //deg.C
    {
      // calculate Halite + Liquid + Vapor coexist composition
      double P_VLH = eos.P_VaporLiquidHaliteCoexist(T);
      double X_VLH = eos.X_VaporHaliteCoexist(T,P_VLH);
      fout<<T<<" "<<X_VLH<<endl;
    }
    fout.close();
  }
  // Fig 8b
  vector<double> arryT={450, 500, 550}; //bar
  for (size_t i = 0; i < arryT.size(); i++)
  {
    string filename="X_HaliteSaturatedVapor_T"+to_string((int)arryT[i])+"C.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double P = 0; P <= 400; P=P+0.5) //bar
    {
      double X_HaliteSaturatedVapor=eos.X_VaporHaliteCoexist(arryT[i],P); //mol fraction
      fout<<P<<" "<<X_HaliteSaturatedVapor<<endl;
    }
    fout.close();
  }
  std::cout<<"Test halite satureated vapor composition end\n"<<endl;
}

void test_Salinity_VaporLiquidCoexist_LiquidBranch()
{
  std::cout<<"Test composition on liquid branch of vapor liquid coexist surface start"<<endl;
  H2ONaCl::cH2ONaCl eos;
  vector<double> arryT={200, 300, 350, 375.5, 380, 400, 500, 600, 800, 1000}; //deg.C
  vector<double> arryPmin={9, 40, 70, 80, 100, 100, 200, 100, 0, 0};
  vector<double> arryPmax={16, 90, 170, 220, 220.4, 270, 560, 900, 1500, 2100};
  vector<double> arryPmax2={16, 90, 170, 229, 236, 290, 600, 950, 2000, 2200};
  vector<double> arryDp={0.05, 0.1, 1, 1, 1, 1, 1, 1, 1, 1};
  double dp_refine = 0.0005;
  for (size_t i = 0; i < arryT.size(); i++)
  {
    string filename="X_VL_LiquidBranch_T"+to_string((int)arryT[i])+"C.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    string filename_vaporBranch="X_VL_VaporBranch_T"+to_string((int)arryT[i])+"C.dat";
    ofstream fout_vaporBranch(filename_vaporBranch);
    if(!fout_vaporBranch)
    {
      cout<<"Open file failed: "<<filename_vaporBranch<<endl;
    }
    for (double P = arryPmin[i]; P <= arryPmax[i]; P=P+arryDp[i]) //bar
    {
      double X_VL_LiquidBranch=eos.X_VaporLiquidCoexistSurface_LiquidBranch(arryT[i],P); //mol fraction
      double X_VL_VaporBranch=eos.X_VaporLiquidCoexistSurface_VaporBranch(arryT[i],P); //mol fraction
      fout<<P<<" "<<X_VL_LiquidBranch<<endl;
      fout_vaporBranch<<P<<" "<<X_VL_VaporBranch<<endl;
    }
    for (double P = arryPmax[i]; P <= arryPmax2[i]; P=P+dp_refine) //bar
    {
      double X_VL_LiquidBranch=eos.X_VaporLiquidCoexistSurface_LiquidBranch(arryT[i],P); //mol fraction
      double X_VL_VaporBranch=eos.X_VaporLiquidCoexistSurface_VaporBranch(arryT[i],P); //mol fraction
      fout<<P<<" "<<X_VL_LiquidBranch<<endl;
      fout_vaporBranch<<P<<" "<<X_VL_VaporBranch<<endl;
    }
    fout.close();
    fout_vaporBranch.close();
  }
  std::cout<<"Test composition on liquid branch of vapor liquid coexist surface end\n"<<endl;
}
void test_V_brine_NaCl_lowThighT()
{
  std::cout<<"Test V_extrapol start"<<endl;
  H2ONaCl::cH2ONaCl eos;
  vector<double> arryP={5, 6, 7, 10, 20,200}; //bar
  vector<double> arryX={0}; //wt.% NaCl
  for (size_t i = 0; i < arryP.size(); i++)
  {
    double P = arryP[i];
    for (size_t j = 0; j < arryX.size(); j++)
    {
      double X_mol = eos.Wt2Mol(arryX[j]/100.0);
      string filename="V_brine_NaCl_lowThighT_P"+to_string((int)P)+"bar_X"+to_string((int)arryX[j])+"wt.dat";
      ofstream fout(filename);
      if(!fout)
      {
        cout<<"Open file failed: "<<filename<<endl;
      }
      for (double T = 0; T <= 1000; T=T+1) //deg.C
      {
        double T_star = eos.T_star_V(T, P, X_mol);
        double V_water = H2O::MolarMass / eos.m_water.Rho(T_star, P);// m3/mol

        double V_extrapol = eos.V_extrapol(T,P,X_mol);

        double T_star_sat = eos.T_star_V(T, P, eos.X_HaliteLiquidus(T,P));
        // cout<<T<<" "<<T_star_sat<<endl;
        double V_water_sat = H2O::MolarMass / eos.m_water.Rho(T_star_sat, P);// m3/mol

        double V_NaCl_liquid=NaCl::MolarMass/eos.m_NaCl.Rho_Liquid(T, P);
        if(T<eos.m_NaCl.T_Melting(P))V_NaCl_liquid=NAN;
        double V_NaCl_solid=NaCl::MolarMass/eos.m_NaCl.Rho_Solid(T, P);
        // if(T<eos.m_NaCl.T_Melting(P))V_NaCl_liquid=NAN;

        fout<<T<<" "<<V_water<<" "<<V_water_sat<<" "<<V_NaCl_liquid<<" "<<T_star_sat<<endl;
      }
      fout.close();
    }
  }
  
  std::cout<<"Test V_extrapol end\n"<<endl;
}
void testT_V_star()
{
  std::cout<<"Test T_V_star start"<<endl;
  H2ONaCl::cH2ONaCl eos;
  vector<double> arryP={1000}; //bar
  vector<double> arryX={0,10,100}; //wt.% NaCl
  for (size_t i = 0; i < arryP.size(); i++)
  {
    double P = arryP[i];
    for (size_t j = 0; j < arryX.size(); j++)
    {
      double X_mol = eos.Wt2Mol(arryX[j]/100.0);
      string filename="V_brine_P"+to_string((int)P)+"bar_X"+to_string((int)arryX[j])+"wt.dat";
      ofstream fout(filename);
      if(!fout)
      {
        cout<<"Open file failed: "<<filename<<endl;
      }
      for (double T = 0; T <= 1000; T=T+1) //deg.C
      {
        double T_star = eos.T_star_V(T, P, X_mol);
        // cout<<T<<" "<<T_star<<endl;
        double V_water = H2O::MolarMass / eos.m_water.Rho(T_star, P);// m3/mol
        // if(X_liquidus>X_VLH)continue;
        fout<<T<<" "<<V_water<<endl;
      }
      fout.close();
    }
  }
  
  std::cout<<"Test T_V_star end\n"<<endl;
}

void test_V_extrapol()
{
  std::cout<<"Test V_extrapol start"<<endl;
  H2ONaCl::cH2ONaCl eos;
  double T=140, P=2, X=0.1;
  double v_extrapol = eos.V_extrapol(T,P,X);
  printf("T: %.1f C, P: %.1f bar, X: %.1f mol frac\n", T, P, X);
  printf("\tH2O boling P: %f bar\n",eos.m_water.BoilingCurve(T));
  printf("\tH2O-NaCl liquidus X: %f mol frac\n",eos.X_HaliteLiquidus(T,P));
  printf("\tV_extrapol: %f\n",v_extrapol);
  std::cout<<"Test V_extrapol end"<<endl;
}

void test_writeCriticalCurve()
{
  eos.writeCriticalCurve();
}
void test_writeHaliteLiquidus()
{
  eos.writeHaliteLiquidusSurface();
}
void test_writePressure_VLH()
{
  eos.writeVaporLiquidHaliteCoexistSurface();
}
void test_NaClMeltingCurve()
{
  eos.writeNaClMeltingCurve();
}
void test_H2OBoilingCurve()
{
  eos.writeH2OBoilingCurve();
}
void test_VaporLiquidHaliteCoexistCurves()
{
  eos.writeVaporLiquidHalite_V_L_H_Curve();
}
void test_VaporLiquidCoexistSurface()
{
  eos.writeVaporLiquidCoexistSurface();
}
int main( int argc, char** argv )
{
  std::cout<<"开始测试计算"<<std::endl;
  // test_CriticalPressure_Composition();
  // test_SublimationBoiling();
  // test_HaliteMelting();
  // test_HaliteLiquidus();
  // test_HaliteSaturatedVaporComposition(); 
  // test_P_VLH();
  // test_Salinity_VaporLiquidCoexist_LiquidBranch();
  // test_water_Curves();
  // test_water_props(130, 200, 2, 10, 0.1, 1, true);
  // testT_V_star();
  // test_V_brine_NaCl_lowThighT();
  // test_V_extrapol();
  // test_NaClH2O_props(0, 800, 1, 1000, 5, 5, true);

  test_writeCriticalCurve();
  test_writeHaliteLiquidus();
  test_writePressure_VLH();
  test_NaClMeltingCurve();
  test_H2OBoilingCurve();
  test_VaporLiquidHaliteCoexistCurves();
  test_VaporLiquidCoexistSurface();


  std::cout<<"测试计算完毕"<<std::endl;
}