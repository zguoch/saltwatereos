#include "H2ONaCl.H"
#include "omp.h"

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

int main( int argc, char** argv )
{
  // test_CriticalPressure_Composition();
  // test_SublimationBoiling();
  // test_HaliteMelting();
  // test_HaliteLiquidus();
  // test_HaliteSaturatedVaporComposition(); 
  // test_P_VLH();
  // test_Salinity_VaporLiquidCoexist_LiquidBranch();
}
void test_P_VLH()
{
  SWEOS::cH2ONaCl eos;
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
  SWEOS::cH2ONaCl eos;
  string filename="HaliteMeltingCurve.dat";
  ofstream fout(filename);
  if(!fout)
  {
    cout<<"Open file failed: "<<filename<<endl;
  }
  for (double P = 0; P <= 5000; P=P+50)
  {
    double T_hm = eos.T_HaliteMelting(P);
    fout<<T_hm<<" "<<P<<endl;
  }
  fout.close();

}
void test_CriticalPressure_Composition()
{
  SWEOS::cH2ONaCl eos;
  // full range
  {
    string filename="HaliteCritical_P_X_fullrange.dat";
    ofstream fout(filename);
    if(!fout)
    {
      cout<<"Open file failed: "<<filename<<endl;
    }
    for (double T = SWEOS::T_Critic_H2O; T <= 1000; T=T+10)
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
    for (double T = SWEOS::T_Critic_H2O; T <= 460; T=T+10)
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
  SWEOS::cH2ONaCl eos;
  string filename="HaliteSublimationBoilingCurve.dat";
  ofstream fout(filename);
  if(!fout)cout<<"Open file failed: "<<filename<<endl;
  
  for (double T = 300; T <= 1100; T=T+10)
  {
    double P_subl = eos.P_HaliteSublimation(T);
    double P_boil = eos.P_HaliteBoiling(T);
    fout<<T<<" "<<P_subl<<" "<<P_boil<<endl;
  }
  fout.close();
}

void test_HaliteLiquidus()
{
  std::cout<<"Test halite liquids start"<<endl;
  SWEOS::cH2ONaCl eos;
  vector<double> arryP={500, 2000, 4000}; //bar
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
  SWEOS::cH2ONaCl eos;
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
  SWEOS::cH2ONaCl eos;
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