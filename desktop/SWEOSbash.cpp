#include "SWEOSbash.h"
// using namespace SWEOSbash;
namespace SWEOSbash
{
  void bash_run(int argc, char** argv)
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
    helpINFO();

    //parse arguments and check 
    cSWEOSarg arg;
    if(!arg.ParseAndCheck(argc, argv)) return;
  }


  cSWEOSarg::cSWEOSarg(/* args */)
  :m_haveD(false), m_haveV(false), m_haveP(false)
  ,m_haveT(false), m_haveX(false), m_haveH(false), m_haveR(false)
  {
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)m_valueR[i][j]=0;
  }

  cSWEOSarg::~cSWEOSarg()
  {
  }

  bool cSWEOSarg::ParseAndCheck(int argc, char** argv)
  {
    int opt; 
    const char *optstring = "D:V:P:T:X:H:R:O:G:t:"; // set argument templete
    while ((opt = getopt(argc, argv, optstring)) != -1) 
    {
      switch (opt)
      {
      case 'D':
        m_haveD=true;
        m_valueD=atoi(optarg);
        if (m_valueD<0 || m_valueD>3)
        {
          cout<<ERROR_COUT<<"option for -D argument must be one of 0, 1, 2, 3. please check -D parameter"<<endl;
          return false;
        }
        break;
      case 'V':
        m_haveV=true;
        m_valueV=optarg;
        for (int i = 0; i < m_valueV.size(); i++)
        {
          if(!(m_valueV[i]=='P' || m_valueV[i]=='T' || m_valueV[i]=='X' || m_valueV[i]=='H'))
          {
            cout<<ERROR_COUT<<"The "<<i+1<<"th variable tag ["<<argv[optind - 1]<<"] cannot be recognized, the supported variables are T, P, X, H"<<endl;
            return false;
          }
        }
        if (m_valueV.size()<1 || m_valueV.size()>3)
        {
          cout<<ERROR_COUT<<"the number of variables cannot exceed three or less than one, it should be one of TPX, THX, T, P, X, H, PT, PX, TX, PH, HX"<<endl;
          return false;
        }
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
      case 'R':
        m_haveR=true;
        m_varueR_str= string_split(optarg,"/");
        if(m_varueR_str.size()%3 != 0)
        {
          cout<<ERROR_COUT<<"range of variables must be a multiple of 3 [min/delta/max], what you set is "<<optarg<<endl;
          return false;
        }
        break;
      default:
        break;
      }
    }
    //check
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
            cout<<ERROR_COUT<<"selected calculation mode is single point PTX, you must specify temperature by -T"<<endl;
            return false;
          }
          if(!m_haveP)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PTX, you must specify pressure by -P"<<endl;
            return false;
          }
          if(!m_haveX)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PTX, you must specify salinity by -X"<<endl;
            return false;
          }
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
          calculateSinglePoint_PTX(m_valueP*1e5, m_valueT, m_valueX);
        }//single point calculation: PHX
        else if(m_valueV=="PHX" || m_valueV=="PXH" || m_valueV=="HPX" || m_valueV=="HXP" || m_valueV=="XPH" || m_valueV=="XPH")
        {
          if(!m_haveH)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PHX, you must specify enthalpy by -H"<<endl;
            return false;
          }
          if(!m_haveP)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PHX, you must specify pressure by -P"<<endl;
            return false;
          }
          if(!m_haveX)
          {
            cout<<ERROR_COUT<<"selected calculation mode is single point PHX, you must specify salinity by -X"<<endl;
            return false;
          }
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
}
