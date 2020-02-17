#include "SWEOSbash.h"
// using namespace SWEOSbash;
namespace SWEOSbash
{
  bool bash_run(int argc, char** argv)
  {
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    if(w.ws_col>119)
    {
        StartText_artASCII();
    }else
    {
        StartText();
    }
    // helpINFO();
    //parse arguments and check 
    cSWEOSarg arg;
    if(!arg.Parse(argc, argv)) return false;
    if(!arg.Validate()) return false;
    return true;
  }


  cSWEOSarg::cSWEOSarg(/* args */)
  :m_haveD(false), m_haveV(false), m_haveP(false)
  ,m_haveT(false), m_haveX(false), m_haveH(false), m_haveR(false),m_haveO(false)
  ,m_valueD(-1), m_valueO(""),m_valueV("")
  {
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)m_valueR[i][j]=0;
  }

  cSWEOSarg::~cSWEOSarg()
  {
  }

  bool cSWEOSarg::Parse(int argc, char** argv)
  {
    if(argc<2)return false; //there is no arguments
    int opt; 
    const char *optstring = "D:V:P:T:X:H:R:O:G:t:vh"; // set argument templete
    int option_index = 0;
    static struct option long_options[] = {
        {"version", no_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, 0, 0}  // to avoid empty input
    };
    int valid_args=0;
    while ((opt = getopt_long(argc, argv, optstring,long_options, &option_index)) != -1) 
    {
      if(opt!='?')valid_args++;
      switch (opt)
      {
      case 'h':
        helpINFO();
        exit(0);
        break;
      case 'v':
        cout<<"Version: "<<VERSION_MAJOR<<"."<<VERSION_MINOR<<endl;
        exit(0);
        break;
      case 'D':
        m_haveD=true;
        m_valueD=atoi(optarg);
        break;
      case 'V':
        m_haveV=true;
        m_valueV=optarg;
        break;
      case 'P':
        m_haveP=true;
        m_valueP=atof(optarg);
        break;
      case 'T':
        m_haveT=true;
        m_valueT=atof(optarg);
        break;
      case 'X':
        m_haveX=true;
        m_valueX=atof(optarg);
        break;
      case 'H':
        m_haveH=true;
        m_valueH=atof(optarg);
        break;
      case 'G':
        m_haveG=true;
        m_valueG=optarg;
        break;
      case 'R':
        m_haveR=true;
        m_varueR_str= string_split(optarg,"/");
        break;
      case 'O':
        m_haveO=true;
        m_valueO=optarg;
      default:
        break;
      }
    }
    if(!(m_haveD && m_haveV))return false;//must have -D and -V arguments
    return true;
  }
  bool cSWEOSarg::Validate()
  {
    //check
    if (m_valueD<0 || m_valueD>3)
    {
      cout<<ERROR_COUT<<"option for -D argument must be one of 0, 1, 2, 3. please check -D parameter"<<endl;
      return false;
    }
    for (int i = 0; i < m_valueV.size(); i++)
    {
      if(!(m_valueV[i]=='P' || m_valueV[i]=='T' || m_valueV[i]=='X' || m_valueV[i]=='H'))
      {
        cout<<ERROR_COUT<<"The option value of -V argument cannot be recognized, the supported variables are T, P, X, H"<<endl;
        return false;
      }
    }
    if (m_valueV.size()<1 || m_valueV.size()>3)
    {
      cout<<ERROR_COUT<<"the number of variables cannot exceed three or less than one, it should be one of TPX, THX, T, P, X, H, PT, PX, TX, PH, HX"<<endl;
      return false;
    }
    if(m_varueR_str.size()%3 != 0)
    {
      cout<<ERROR_COUT<<"range of variables must be a multiple of 3 [min/delta/max]"<<endl;
      return false;
    }
    if(!m_haveD)
    {
      cout<<ERROR_COUT<<"You have to specify the -D argument, the parameter should be one of 0, 1, 2, 3"<<endl;
      return false;
    }
    if(!m_haveV)
    {
      cout<<ERROR_COUT<<"You have to specify the -V argument, the parameter should be one of PTX, PHX, T, P, X, H, PT, PX, TX, PH, HX"<<endl;
      return false;
    }
    switch (m_valueD)
    {
    case CALCULATION_MODE_SINGLEPOINT:
      {
        if(m_valueV.size()!=3)
        {
          cout<<ERROR_COUT<<"if -D set as -D0, the -V support [P,T,X] or [P,H,X]: "<<endl;
          return false;
        }
        //single point calculation: PTX
        if(m_valueV=="PTX" || m_valueV=="PXT" || m_valueV=="TPX" || m_valueV=="TXP" || m_valueV=="XPT" || m_valueV=="XPT")
        {
          if(!m_haveT)
          {
            cout<<WARN_COUT<<"selected calculation mode is single point PTX, you must specify temperature by -T"<<endl;
          }
          if(!m_haveP)
          {
            cout<<WARN_COUT<<"selected calculation mode is single point PTX, you must specify pressure by -P"<<endl;
          }
          if(!m_haveX)
          {
            cout<<WARN_COUT<<"selected calculation mode is single point PTX, you must specify salinity by -X"<<endl;
          }
          // determin single point calculation, multi-points calculation or exit
          if(m_haveT && m_haveP && m_haveX)//single point calculation
          {
            if(m_haveG)cout<<WARN_COUT<<"have -T, -P, -X option values, ignore -G argument"<<endl;
            //check range
            if(m_valueP>SWEOS::PMAX/1e5 || m_valueP<SWEOS::PMIN/1e5)
            {
              cout<<ERROR_COUT<<"-P specify pressure ="<<m_valueP<<"bar = "<<m_valueP*1e5<<" Pa, out of range ["<<SWEOS::PMIN<<", "<<SWEOS::PMAX<<"] Pa"<<endl;
              return false;
            }
            if(m_valueT>SWEOS::TMAX-SWEOS::Kelvin || m_valueT<SWEOS::TMIN-SWEOS::Kelvin)
            {
              cout<<ERROR_COUT<<"-T specify temperature ="<<m_valueT<<"°C = "<<m_valueT+SWEOS::Kelvin<<" K, out of range ["<<SWEOS::TMIN<<", "<<SWEOS::TMAX<<"] K"<<endl;
              return false;
            }
            if(m_valueX>SWEOS::XMAX || m_valueX<SWEOS::XMIN)
            {
              cout<<ERROR_COUT<<"-X specify salinity ="<<m_valueX<<" wt. % NaCl, out of range ["<<SWEOS::XMIN<<", "<<SWEOS::XMAX<<"] wt. % NaCl"<<endl;
              return false;
            }
            calculateSinglePoint_PTX(m_valueP*1e5, m_valueT+SWEOS::Kelvin, m_valueX);
          }else if(m_haveG)//not T, P, X value, but have intput file name
          {
            cout<<WARN_COUT<<"You specify a input file for multi-points calculation\n"
            <<"Please make sure your input file with 3 columns in order of "<<m_valueV<<endl;
            calculateMultiPoints_PTX_PHX(m_valueV,m_valueG, m_valueO,"T");
          }else
          {
            cout<<ERROR_COUT<<"There neither full -T, -P, -X options nor -G argument, swEOS will exit"<<endl;
            exit(0);
          }
        }//single point calculation: PHX
        else if(m_valueV=="PHX" || m_valueV=="PXH" || m_valueV=="HPX" || m_valueV=="HXP" || m_valueV=="XPH" || m_valueV=="XPH")
        {
          if(!m_haveH)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PHX, you must specify enthalpy by -H"<<endl;
          }
          if(!m_haveP)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PHX, you must specify pressure by -P"<<endl;
          }
          if(!m_haveX)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PHX, you must specify salinity by -X"<<endl;
          }
          
          // determin single point calculation, multi-points calculation or exit
          if(m_haveH && m_haveP && m_haveX)//single point calculation
          {
            //check range
            if(m_valueP>SWEOS::PMAX/1e5 || m_valueP<SWEOS::PMIN/1e5)
            {
              cout<<ERROR_COUT<<"-P specify pressure ="<<m_valueP<<"bar = "<<m_valueP*1e5<<" Pa, out of range ["<<SWEOS::PMIN<<", "<<SWEOS::PMAX<<"] Pa"<<endl;
              return false;
            }
            if(m_valueX>SWEOS::XMAX || m_valueX<SWEOS::XMIN)
            {
              cout<<ERROR_COUT<<"-X specify salinity ="<<m_valueX<<" wt. % NaCl, out of range ["<<SWEOS::XMIN<<", "<<SWEOS::XMAX<<"] wt. % NaCl"<<endl;
              return false;
            }
            calculateSinglePoint_PHX(m_valueP*1e5, m_valueH*1000, m_valueX);
          }else if(m_haveG)//not T, P, X value, but have intput file name
          {
            cout<<WARN_COUT<<"You specify a input file for multi-points calculation\n"
            <<"Please make sure your input file with 3 columns in order of "<<m_valueV<<endl;
            calculateMultiPoints_PTX_PHX(m_valueV,m_valueG, m_valueO,"H");
          }else
          {
            cout<<ERROR_COUT<<"There neither full -H, -P, -X options nor -G argument, swEOS will exit"<<endl;
            exit(0);
          }
        }else
        {
          cout<<ERROR_COUT<<"unknown -D parameter: "<<m_valueV<<endl;
          return false;
        }
        
      }
      break;
    
    default:
      break;
    }
    return true;
  }
  vector<string> cSWEOSarg::string_split (string s, string delimiter) {
      size_t pos_start = 0, pos_end, delim_len = delimiter.length();
      string token;
      vector<string> res;
      while ((pos_end = s.find (delimiter, pos_start)) != string::npos) 
      {
          token = s.substr (pos_start, pos_end - pos_start);
          pos_start = pos_end + delim_len;
          res.push_back (token);
      }
      res.push_back (s.substr (pos_start));
      return res;
  }
  SWEOS::PROP_H2ONaCl calculateSinglePoint_PTX(double P, double T, double X, bool isCout)
  {
    SWEOS::cH2ONaCl eos;
    eos.prop_pTX(P, T+SWEOS::Kelvin, X);
    if(isCout) 
    {
      eos.setColorPrint(true);
      cout<<eos<<endl;
    }
    return eos.m_prop;
  }
  SWEOS::PROP_H2ONaCl calculateSinglePoint_PHX(double P, double H, double X, bool isCout)
  {
    SWEOS::cH2ONaCl eos;
    eos.prop_pHX(P, H, X);
    if(isCout)
    {
      eos.setColorPrint(true);
      cout<<eos<<endl;
    }
    return eos.m_prop;
  }
  bool calculateMultiPoints_PTX_PHX(string valueV, string filePTX, string outFile, string isT_H)
  {
    ifstream fin(filePTX);
    if(!fin)
    {
      cout<<ERROR_COUT<<"Open file failed, please check -G argument, the file name specified by -G is "<<COLOR_RED<<filePTX<<COLOR_DEFAULT<<endl;
      exit(0);
    }
    vector<double>P, T_H, X;
    double p,t,x;
    if(valueV=="P"+isT_H+"X")
    {
      while (!fin.eof())
      {
        fin>>p>>t>>x;
        P.push_back(p);
        T_H.push_back(t);
        X.push_back(x);
      }
    }else if(valueV=="PX"+isT_H+"")
    {
      while (!fin.eof())
      {
        fin>>p>>x>>t;
        P.push_back(p);
        T_H.push_back(t);
        X.push_back(x);
      }
    }else if(valueV==""+isT_H+"PX")
    {
      while (!fin.eof())
      {
        fin>>t>>p>>x;
        P.push_back(p);
        T_H.push_back(t);
        X.push_back(x);
      }
    }else if(valueV==""+isT_H+"XP")
    {
      while (!fin.eof())
      {
        fin>>t>>x>>p;
        P.push_back(p);
        T_H.push_back(t);
        X.push_back(x);
      }
    }else if(valueV=="XP"+isT_H+"")
    {
      while (!fin.eof())
      {
        fin>>x>>p>>t;
        P.push_back(p);
        T_H.push_back(t);
        X.push_back(x);
      }
    }else if(valueV=="X"+isT_H+"P")
    {
      while (!fin.eof())
      {
        fin>>x>>t>>p;
        P.push_back(p);
        T_H.push_back(t);
        X.push_back(x);
      }
    }else
    {
      cout<<ERROR_COUT<<"The -V argument must be one of P"+isT_H+"X, PX"+isT_H+
                        ", "+isT_H+"PX, "+isT_H+"XP, XP"+isT_H+", X"+isT_H+"P when -D0\n"
      <<"The -V option you set is "<<COLOR_RED<<valueV<<COLOR_DEFAULT<<" which is not supported"<<endl;
      exit(0);
    }
    fin.close();
    //cout or fout
    ofstream fout(outFile);
    SWEOS::cH2ONaCl eos;
    if(fout)
    {
      fout<<"T(C), P(bar), X, Phase Index, Phase Region, "
          <<"Bulk density(kg/m3), Liquid density(kg/m3), Vapour density(kg/m3), Halite density(kg/m3), "
          <<"Bulk enthalpy(kJ/kg), Liquid enthalpy(kJ/kg), Vapour enthalpy(kJ/kg), Halite enthalpy(kJ/kg), "
          <<"Liquid Saturation, Vapour Saturation, Halite Saturation, "
          <<"Liquid viscosity, Vapour viscosity, "
          <<"Liquid salinity, Vapour salinity"
          <<endl;
      for(int i=0;i<T_H.size();i++)
      {
        if(isT_H=="T")
        {
          eos.prop_pTX(P[i]*1e5, T_H[i]+SWEOS::Kelvin, X[i]);
        }else if(isT_H=="H")
        {
          eos.prop_pHX(P[i]*1e5, T_H[i]*1000, X[i]);
        }
        fout<<eos.m_prop.T+SWEOS::Kelvin<<", "<<P[i]<<", "<<X[i]<<", "<<eos.m_prop.Region<<", "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<", "
            <<eos.m_prop.Rho<<", "<<eos.m_prop.Rho_l<<", "<<eos.m_prop.Rho_v<<", "<<eos.m_prop.Rho_h<<", "
            <<eos.m_prop.H<<", "<<eos.m_prop.H_l<<", "<<eos.m_prop.H_v<<", "<<eos.m_prop.H_h<<", "
            <<eos.m_prop.S_l<<", "<<eos.m_prop.S_v<<", "<<eos.m_prop.S_h<<", "
            <<eos.m_prop.Mu_l<<", "<<eos.m_prop.Mu_v<<", "
            <<eos.m_prop.X_l<<", "<<eos.m_prop.X_v
            <<endl;
      }
      fout.close();
    }else
    {
      if(!fout)cout<<WARN_COUT<<"The output file name is specified by -O argument, but open file failed: "<<outFile
                   <<"\nPrint in terminal instead"<<endl;
      cout<<"T(C), P(bar), X, Phase Index, Phase Region, "
          <<"Bulk density(kg/m3), Liquid density(kg/m3), Vapour density(kg/m3), Halite density(kg/m3), "
          <<"Bulk enthalpy(kJ/kg), Liquid enthalpy(kJ/kg), Vapour enthalpy(kJ/kg), Halite enthalpy(kJ/kg), "
          <<"Liquid Saturation, Vapour Saturation, Halite Saturation, "
          <<"Liquid viscosity, Vapour viscosity, "
          <<"Liquid salinity, Vapour salinity"
          <<endl;
      for(int i=0;i<T_H.size();i++)
      {
        if(isT_H=="T")
        {
          eos.prop_pTX(P[i]*1e5, T_H[i]+SWEOS::Kelvin, X[i]);
        }else if(isT_H=="H")
        {
          eos.prop_pHX(P[i]*1e5, T_H[i]*1000, X[i]);
        }
        cout<<eos.m_prop.T+SWEOS::Kelvin<<", "<<P[i]<<", "<<X[i]<<", "<<eos.m_prop.Region<<", "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<", "
            <<eos.m_prop.Rho<<", "<<eos.m_prop.Rho_l<<", "<<eos.m_prop.Rho_v<<", "<<eos.m_prop.Rho_h<<", "
            <<eos.m_prop.H<<", "<<eos.m_prop.H_l<<", "<<eos.m_prop.H_v<<", "<<eos.m_prop.H_h<<", "
            <<eos.m_prop.S_l<<", "<<eos.m_prop.S_v<<", "<<eos.m_prop.S_h<<", "
            <<eos.m_prop.Mu_l<<", "<<eos.m_prop.Mu_v<<", "
            <<eos.m_prop.X_l<<", "<<eos.m_prop.X_v
            <<endl;
      }
    }
    return true;
  }
}
