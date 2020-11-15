#include "H2ONaCl.H"
#include "omp.h"

// 1. Halite liquidus
void test_HaliteLiquidus();

int main( int argc, char** argv )
{
  test_HaliteLiquidus();
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