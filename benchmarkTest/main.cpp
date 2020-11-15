#include "H2ONaCl.H"
#include "omp.h"

// 1. Halite melting curve, Fig. 3 of Driesner and Heinrich(2007)
void test_HaliteMelting();
// 2. Vapor pressure curves (sublimation and boling curves), Fig. 4 of Driesner and Heinrich(2007)
void test_SublimationBoiling();
// 3. Critical pressure and salinity , Fig. 5 of Driesner and Heinrich(2007)
void test_CriticalPressure_Composition();
// 1. Halite liquidus, Fig. 7 of Driesner and Heinrich(2007)
void test_HaliteLiquidus();

int main( int argc, char** argv )
{
  test_CriticalPressure_Composition();
  // test_SublimationBoiling();
  // test_HaliteMelting();
  // test_HaliteLiquidus();
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