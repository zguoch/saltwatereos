#include "SWEOSbash.h"
// using namespace SWEOSbash;
namespace SWEOSbash
{
  bool bash_run(int argc, char** argv)
  {
    #ifdef _WIN32
        
    #else
      StartText();
      // struct winsize w;
      // ioctl(0, TIOCGWINSZ, &w);
      // if(w.ws_col>119)
      // {
      //     StartText_artASCII();
      // }else
      // {
      //     StartText();
      // }
    #endif
    
    helpINFO();
    //parse arguments and check 
    cSWEOSarg arg;
    if(!arg.Parse(argc, argv)) return false;
    if(!arg.Validate()) return false;
    return true;
  }


  cSWEOSarg::cSWEOSarg(/* args */)
  :m_haveD(false), m_haveV(false), m_haveP(false)
  ,m_haveT(false), m_havet(false), m_haveX(false), m_haveH(false), m_haveR(false),m_haveO(false)
  ,m_valueD(-1),m_threadNumOMP(omp_get_max_threads()), m_valueO(""),m_valueV("")
  {
    for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)m_valueR[i][j]=0;
  }

  cSWEOSarg::~cSWEOSarg()
  {
  }
  bool isNum(string str)
  {
      stringstream sin(str);
      double d;
      char c;
      if(!(sin >> d))
          return false;
      if (sin >> c)
          return false;
      return true;
  }
  bool cSWEOSarg::GetOptionValue(int opt, char* optarg, double& value)
  {
    string optarg_str=optarg;
    if(isNum(optarg_str))
    {
      value=atof(optarg);
    }else
    { 
      char optCh=opt;
      cout<<ERROR_COUT<<"Option of -"<<optCh<<" argument is empty or cannot be recognized"<<endl;
      return false;
    }
    return true;
  }
  
  #ifdef _WIN32
    bool cSWEOSarg::Parse(int argc, char** argv)
    {
      return true;
    }
  #else
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
      double doubleOptValue;
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
          if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
          m_valueD=(int)doubleOptValue;
          break;
        case 't':
          m_havet=true;
          if(!GetOptionValue(opt, optarg, doubleOptValue))return false;
          m_threadNumOMP=(int)doubleOptValue;
          if(m_threadNumOMP>omp_get_max_threads())m_threadNumOMP=omp_get_max_threads();
          if(m_threadNumOMP<1)m_threadNumOMP=1;
          break;
        case 'V':
          m_haveV=true;
          m_valueV=optarg;
          break;
        case 'P':
          m_haveP=true;
          if(!GetOptionValue(opt, optarg, m_valueP))return false;
          break;
        case 'T':
          m_haveT=true;
          if(!GetOptionValue(opt, optarg, m_valueT))return false;
          break;
        case 'X':
          m_haveX=true;
          if(!GetOptionValue(opt, optarg, m_valueX))return false;
          break;
        case 'H':
          m_haveH=true;
          if(!GetOptionValue(opt, optarg, m_valueH))return false;
          break;
        case 'G':
          m_haveG=true;
          m_valueG=optarg;
          break;
        case 'R':
          m_haveR=true;
          m_valueR_str= string_split(optarg,"/");
          break;
        case 'O':
          m_haveO=true;
          m_valueO=optarg;
          break;
        default:
          break;
        }
      }
      if(!(m_haveD && m_haveV))return false;//must have -D and -V arguments
      return true;
    }
  #endif

  bool cSWEOSarg::Validate()
  {
    //check required arguments
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
    if(m_valueR_str.size()%3 != 0 && m_valueR_str.size()>9)
    {
      cout<<ERROR_COUT<<"Option of -R argument must be a multiple of 3 and <=9, in format of [min/delta/max]"<<endl;
      return false;
    }else //set m_valueR
    {
      for (size_t i = 0; i < m_valueR_str.size(); i++)
      {
        if(isNum(m_valueR_str[i]))
        {
          m_valueR[i/3][i%3]=atof(m_valueR_str[i].c_str());
        }else
        {
          cout<<ERROR_COUT<<"The "<<i+1<<"th value in -R option is not a number: "<<COLOR_RED<<m_valueR_str[i]<<COLOR_DEFAULT<<endl;
          return false;
        }
      }
      
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
        return Validate_0D();
      }
      break;
    case CALCULATION_MODE_ONEDIMENSION:
      {
        return Validate_1D();
      }
      break;
    case CALCULATION_MODE_TWODIMENSION:
      {
        return Validate_2D();
      }
      break;
    case CALCULATION_MODE_THREEDIMENSION:
      {
        return Validate_3D();
      }
      break;
    default:
      break;
    }
    return true;
  }
  bool cSWEOSarg::Validate_2D()
  {
    if(m_valueV.size()!=2)
    {
      cout<<ERROR_COUT<<"if -D set as -D2, the -V option must be one of -VPT, -VPX, -VTX, -VPH, -VHX"<<endl;
      return false;
    }
    if (!m_haveR)
    {
      cout<<ERROR_COUT<<"You set -D2 and -V"<<m_valueV<<", then you must set -R for range of "
          <<COLOR_GREEN<<m_valueV[0]<<COLOR_DEFAULT<<" and "<<COLOR_GREEN<<m_valueV[1]<<COLOR_DEFAULT<<" in format of -R"
          <<COLOR_GREEN<<m_valueV[0]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[0]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[0]<<"max/"
          <<COLOR_GREEN<<m_valueV[1]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[1]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[1]<<"max"<<COLOR_DEFAULT
          <<endl;
      return false;
    }
    if(m_valueR_str.size()!=6)
    {
      cout<<ERROR_COUT<<"You set -D2 and -V"<<m_valueV<<", then you must set -R for range of "
          <<COLOR_GREEN<<m_valueV[0]<<COLOR_DEFAULT<<" and "<<COLOR_GREEN<<m_valueV[1]<<COLOR_DEFAULT<<" in format of -R"
          <<COLOR_GREEN<<m_valueV[0]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[0]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[0]<<"max/"
          <<COLOR_GREEN<<m_valueV[1]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[1]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[1]<<"max"<<COLOR_DEFAULT
          <<endl;
      cout<<ERROR_COUT<<"Option of -R must be 6 values, but what you set is "<<COLOR_RED;
      for(int i=0;i<m_valueR_str.size();i++)cout<<m_valueR_str[i]<<" ";
      cout<<COLOR_DEFAULT<<endl;

      return false;
    }
    if(!m_haveO || m_valueO=="")
    {
      cout<<WARN_COUT<<"You forget to set output file name through -O argument, but doesn't matter, it is reseted as "
          <<m_valueV<<".vtk"<<endl;
      m_valueO=m_valueV+".vtk";
    }
    if(m_valueV=="TP" || m_valueV=="PT")// fixed X
    {
      if(!(m_haveX && CheckRange_X(m_valueX)))
      {
        cout<<ERROR_COUT<<"Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument"<<endl;
        return false;
      }
      int ind_T=0, ind_P=1;
      if(m_valueV=="PT")
      {
        ind_P=0; ind_T=1;
      }
      double rangeT[2]={m_valueR[ind_T][0], m_valueR[ind_T][2]};
      double rangeP[2]={m_valueR[ind_P][0], m_valueR[ind_P][2]};
      if(!CheckRanges_T(rangeT)) return false;
      if(!CheckRanges_P(rangeP)) return false;
      //calculate
      vector<double> arrT= linspace(m_valueR[ind_T][0], m_valueR[ind_T][2], m_valueR[ind_T][1]);
      vector<double> arrP= linspace(m_valueR[ind_P][0], m_valueR[ind_P][2], m_valueR[ind_P][1]);
      vector<double> arrX; arrX.push_back(m_valueX);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrT.size()*arrP.size());
      
      MultiProgressBar multibar(arrP.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"2D calculation using "<<m_threadNumOMP<<" threads, T ∈ ["
          <<rangeT[0]<<", "<<rangeT[1]<<"] °C, P ∈ ["
          <<rangeP[0]<<", "<<rangeP[1]<<"] bar, fixed salinity X="
          <<m_valueX<<" "
          <<"\n"<<endl;
      int lenT = (int)(arrT.size());
      int lenP = (int)(arrP.size());
      #pragma omp parallel for shared(arrT, arrP, arrX, props, lenT)
      for (int j = 0; j < lenP; j++)
      {
        for (int k = 0; k < lenT; k++)
        {
          SWEOS::cH2ONaCl eos;
          eos.prop_pTX(arrP[j]*1e5, arrT[k]+SWEOS::Kelvin, arrX[0]);
          props[k+j*lenT]=eos.m_prop;
        }
        #pragma omp critical
        multibar.Update();
      } 
      Write2D3DResult(arrT, arrP, arrX, props, m_valueO, "Temperature (°C)", "Pressure (bar)", "Salinity");
        
    }else if(m_valueV=="PX" || m_valueV=="XP")
    {
      if(!(m_haveT && CheckRange_T(m_valueT)))
      {
        cout<<ERROR_COUT<<"Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed temperature value by -T argument"<<endl;
        return false;
      }
      int ind_X=0, ind_P=1;
      if(m_valueV=="PX")
      {
        ind_P=0; ind_X=1;
      }
      double rangeX[2]={m_valueR[ind_X][0], m_valueR[ind_X][2]};
      double rangeP[2]={m_valueR[ind_P][0], m_valueR[ind_P][2]};
      if(!CheckRanges_X(rangeX)) return false;
      if(!CheckRanges_P(rangeP)) return false;
      //calculate
      vector<double> arrX= linspace(m_valueR[ind_X][0], m_valueR[ind_X][2], m_valueR[ind_X][1]);
      vector<double> arrP= linspace(m_valueR[ind_P][0], m_valueR[ind_P][2], m_valueR[ind_P][1]);
      vector<double> arrT; arrT.push_back(m_valueT);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrP.size()*arrX.size());
      MultiProgressBar multibar(arrP.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"2D calculation using "<<m_threadNumOMP<<" threads, X ∈ ["
          <<rangeX[0]<<", "<<rangeX[1]<<"] , P ∈ ["
          <<rangeP[0]<<", "<<rangeP[1]<<"] bar, fixed temperature T="
          <<m_valueT<<" °C "
          <<"\n"<<endl;
      int lenX=arrX.size();
      int lenP = (int)(arrP.size());
      #pragma omp parallel for shared(arrT, arrP, arrX, props, lenX)
      for (int j = 0; j < lenP; j++)
      {
        for (int k = 0; k < lenX; k++)
        {
          SWEOS::cH2ONaCl eos;
          eos.prop_pTX(arrP[j]*1e5, arrT[0]+SWEOS::Kelvin, arrX[k]);
          props[k+j*lenX]=eos.m_prop;
        }
        #pragma omp critical
        multibar.Update();
      } 
      Write2D3DResult(arrX, arrP, arrT, props, m_valueO, "Salinity", "Pressure (bar)", "Temperature (°C)");
    }else if(m_valueV=="TX" || m_valueV=="XT")
    {
      if(!(m_haveP && CheckRange_P(m_valueP)))
      {
        cout<<ERROR_COUT<<"Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument"<<endl;
        return false;
      }
      int ind_X=0, ind_T=1;
      if(m_valueV=="TX")
      {
        ind_T=0; ind_X=1;
      }
      double rangeX[2]={m_valueR[ind_X][0], m_valueR[ind_X][2]};
      double rangeT[2]={m_valueR[ind_T][0], m_valueR[ind_T][2]};
      if(!CheckRanges_X(rangeX)) return false;
      if(!CheckRanges_T(rangeT)) return false;
      //calculate
      vector<double> arrX= linspace(m_valueR[ind_X][0], m_valueR[ind_X][2], m_valueR[ind_X][1]);
      vector<double> arrT= linspace(m_valueR[ind_T][0], m_valueR[ind_T][2], m_valueR[ind_T][1]);
      vector<double> arrP; arrP.push_back(m_valueP);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrT.size()*arrX.size());
      MultiProgressBar multibar(arrX.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"2D calculation using "<<m_threadNumOMP<<" threads, X ∈ ["
          <<rangeX[0]<<", "<<rangeX[1]<<"] , T ∈ ["
          <<rangeT[0]<<", "<<rangeT[1]<<"] °C, fixed pressure P="
          <<m_valueP<<" bar "
          <<"\n"<<endl;
      int lenT=arrT.size();
      int lenX = (int)(arrX.size());
      #pragma omp parallel for shared(arrT, arrP, arrX, props, lenT)
      for (int j = 0; j < lenX; j++)
      {
        for (int k = 0; k < lenT; k++)
        {
          SWEOS::cH2ONaCl eos;
          eos.prop_pTX(arrP[0]*1e5, arrT[k]+SWEOS::Kelvin, arrX[j]);
          props[k+j*lenT]=eos.m_prop;
        }
        #pragma omp critical
        multibar.Update();
      } 
      Write2D3DResult(arrT, arrX, arrP, props, m_valueO, "Temperature (°C)", "Salinity", "Pressure (bar)");
    }else if(m_valueV=="PH" || m_valueV=="HP")
    {
      if(!(m_haveX && CheckRange_X(m_valueX)))
      {
        cout<<ERROR_COUT<<"Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument"<<endl;
        return false;
      }
      int ind_H=0, ind_P=1;
      if(m_valueV=="PH")
      {
        ind_P=0; ind_H=1;
      }
      double rangeH[2]={m_valueR[ind_H][0], m_valueR[ind_H][2]};
      double rangeP[2]={m_valueR[ind_P][0], m_valueR[ind_P][2]};
      if(!CheckRanges_H_P(rangeH[0], rangeH[1], rangeP, m_valueX)) return false;
      if(!CheckRanges_P(rangeP)) return false;
      //calculate
      vector<double> arrH= linspace(m_valueR[ind_H][0], m_valueR[ind_H][2], m_valueR[ind_H][1]);
      vector<double> arrP= linspace(m_valueR[ind_P][0], m_valueR[ind_P][2], m_valueR[ind_P][1]);
      vector<double> arrX; arrX.push_back(m_valueX);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrH.size()*arrP.size());
      
      MultiProgressBar multibar(arrP.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"2D calculation using "<<m_threadNumOMP<<" threads, H ∈ ["
          <<rangeH[0]<<", "<<rangeH[1]<<"] kJ/kg, P ∈ ["
          <<rangeP[0]<<", "<<rangeP[1]<<"] bar, fixed salinity X="
          <<m_valueX<<" "
          <<"\n"<<endl;
      int lenH=arrH.size();
      int lenP = (int)(arrP.size());
      #pragma omp parallel for shared(arrH, arrP, arrX, props, lenH)
      for (int j = 0; j < lenP; j++)
      {
        for (int k = 0; k < lenH; k++)
        {
          SWEOS::cH2ONaCl eos;
          eos.prop_pHX(arrP[j]*1e5, arrH[k]*1000.0, arrX[0]);
          props[k+j*lenH]=eos.m_prop;
        }
        #pragma omp critical
        multibar.Update();
      } 
      Write2D3DResult(arrH, arrP, arrX, props, m_valueO, "Enthalpy (kJ/kg)", "Pressure (bar)", "Salinity");
    }else if(m_valueV=="HX" || m_valueV=="XH")
    {
      if(!(m_haveP && CheckRange_P(m_valueP)))
      {
        cout<<ERROR_COUT<<"Selected calculation mode is 2D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument"<<endl;
        return false;
      }
      int ind_H=0, ind_X=1;
      if(m_valueV=="XH")
      {
        ind_X=0; ind_H=1;
      }
      double rangeH[2]={m_valueR[ind_H][0], m_valueR[ind_H][2]};
      double rangeX[2]={m_valueR[ind_X][0], m_valueR[ind_X][2]};
      if(!CheckRanges_H_X(rangeH[0], rangeH[1], rangeX, m_valueP)) return false;
      if(!CheckRanges_X(rangeX)) return false;
      //calculate
      vector<double> arrH= linspace(m_valueR[ind_H][0], m_valueR[ind_H][2], m_valueR[ind_H][1]);
      vector<double> arrX= linspace(m_valueR[ind_X][0], m_valueR[ind_X][2], m_valueR[ind_X][1]);
      vector<double> arrP; arrP.push_back(m_valueP);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrH.size()*arrX.size());
      
      MultiProgressBar multibar(arrX.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"2D calculation using "<<m_threadNumOMP<<" threads, H ∈ ["
          <<rangeH[0]<<", "<<rangeH[1]<<"] kJ/kg, X ∈ ["
          <<rangeX[0]<<", "<<rangeX[1]<<"], fixed pressure P="
          <<m_valueP<<" bar"
          <<"\n"<<endl;
      int lenH=arrH.size();
      int lenX = (int)(arrX.size());
      #pragma omp parallel for shared(arrH, arrP, arrX, props, lenH)
      for (int j = 0; j < lenX; j++)
      {
        for (int k = 0; k < lenH; k++)
        {
          SWEOS::cH2ONaCl eos;
          eos.prop_pHX(arrP[0]*1e5, arrH[k]*1000.0, arrX[j]);
          props[k+j*lenH]=eos.m_prop;
        }
        #pragma omp critical
        multibar.Update();
      } 
      Write2D3DResult(arrH, arrX, arrP, props, m_valueO, "Enthalpy (kJ/kg)", "Salinity", "Pressure (bar)");
    }else
    {
      cout<<ERROR_COUT<<"Unrecognized -V parameter for two-dimension calculation: -V"<<m_valueV<<endl;
      cout<<"\tAvailable options are: -VPT, -VPX, -VTX, -VPH, -VHX"<<endl;
      return false;
    }
    return true;
  }
  bool Write2D3DResult(vector<double> x, vector<double> y, vector<double> z, vector<SWEOS::PROP_H2ONaCl> props, 
                       string outFile, string xTitle, string yTitle, string zTitle, bool isWritePy)
  {
    string extname;
    string fname_pyScript;
    vector<string> tmp=string_split(outFile,".");
    if(tmp.size()>=1)
    {
      extname=tmp[tmp.size()-1];
    }
    if(extname=="vtk")
    {
      SWEOS::cH2ONaCl eos;
      eos.writeProps2VTK(x,y,z,props, outFile, isWritePy, xTitle, yTitle, zTitle);
    }else
    { 
      cout<<WARN_COUT<<"Unrecognized format: "<<outFile<<endl;
      cout<<COLOR_GREEN<<"Write results into vtk file format"<<COLOR_DEFAULT<<endl;
      string newfilename="";
      for (size_t i = 0; i < tmp.size()-1; i++)
      {
        newfilename+=tmp[i];
      }
      outFile=newfilename+".vtk";
      SWEOS::cH2ONaCl eos;
      eos.writeProps2VTK(x,y,z,props, outFile, isWritePy, xTitle, yTitle, zTitle);
    }
    cout<<COLOR_BLUE<<"Results have been saved to file: "<<outFile<<endl;
    cout<<COLOR_BLUE<<"Paraview-python script is generated as : "<<outFile+".py"<<endl;
    string cmd_pv=" paraview --script="+outFile+".py";
    cout<<"You can use command of "<<COLOR_GREEN<<cmd_pv<<COLOR_DEFAULT<<" to present result in paraview"<<endl;
    return true;
  }
  bool cSWEOSarg::Validate_3D()
  {
    if(m_valueV.size()!=3)
    {
      cout<<ERROR_COUT<<"if -D set as -D3, the -V option must be one of -VPTX, -VPHX"<<endl;
      return false;
    }
    if (!m_haveR)
    {
      cout<<ERROR_COUT<<"You set -D3 and -V"<<m_valueV<<", then you must set -R for range of "
          <<COLOR_GREEN<<m_valueV[0]<<COLOR_DEFAULT<<", "<<COLOR_GREEN<<m_valueV[1]<<COLOR_DEFAULT<<" and "<<COLOR_GREEN<<m_valueV[2]<<COLOR_DEFAULT<<" in format of -R"
          <<COLOR_GREEN<<m_valueV[0]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[0]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[0]<<"max/"
          <<COLOR_GREEN<<m_valueV[1]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[1]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[1]<<"max"<<COLOR_DEFAULT
          <<COLOR_GREEN<<m_valueV[2]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[2]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[1]<<"max"<<COLOR_DEFAULT
          <<endl;
      return false;
    }
    if(m_valueR_str.size()!=9)
    {
      cout<<ERROR_COUT<<"You set -D3 and -V"<<m_valueV<<", then you must set -R for range of "
          <<COLOR_GREEN<<m_valueV[0]<<COLOR_DEFAULT<<", "<<COLOR_GREEN<<m_valueV[1]<<COLOR_DEFAULT<<" and "<<COLOR_GREEN<<m_valueV[2]<<COLOR_DEFAULT<<" in format of -R"
          <<COLOR_GREEN<<m_valueV[0]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[0]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[0]<<"max/"
          <<COLOR_GREEN<<m_valueV[1]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[1]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[1]<<"max"<<COLOR_DEFAULT
          <<COLOR_GREEN<<m_valueV[2]<<"min"<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<"d"<<m_valueV[2]<<COLOR_DEFAULT<<"/"<<COLOR_GREEN<<m_valueV[1]<<"max"<<COLOR_DEFAULT
          <<endl;
      cout<<ERROR_COUT<<"Option of -R must be 9 values, but what you set is "<<COLOR_RED;
      for(int i=0;i<m_valueR_str.size();i++)cout<<m_valueR_str[i]<<" ";
      cout<<COLOR_DEFAULT<<endl;

      return false;
    }
    if(!m_haveO || m_valueO=="")
    {
      cout<<WARN_COUT<<"You forget to set output file name through -O argument, but doesn't matter, it is reseted as "
          <<m_valueV<<".vtk"<<endl;
      m_valueO=m_valueV+".vtk";
    }
    if(m_valueV=="PTX" || m_valueV=="PXT" || m_valueV=="TPX" || m_valueV=="TXP" || m_valueV=="XPT" || m_valueV=="XTP")
    {
      int indP=0, indT=1, indX=2;
      if(m_valueV=="PXT")
      {
        indP=0; indX=1; indT=2;
      }else if(m_valueV=="TPX")
      {
        indT=0; indP=1; indX=2;
      }else if(m_valueV=="TXP")
      {
        indT=0; indX=1; indP=2;
      }else if(m_valueV=="XPT")
      {
        indX=0; indP=1; indT=2;
      }else if(m_valueV=="XTP")
      {
        indX=0; indT=1; indP=2;
      }
      double rangeT[2]={m_valueR[indT][0], m_valueR[indT][2]};
      double rangeP[2]={m_valueR[indP][0], m_valueR[indP][2]};
      double rangeX[2]={m_valueR[indX][0], m_valueR[indX][2]};
      if(!CheckRanges_T(rangeT)) return false;
      if(!CheckRanges_P(rangeP)) return false;
      if(!CheckRanges_X(rangeX)) return false;
      //calculate
      vector<double> arrT= linspace(m_valueR[indT][0], m_valueR[indT][2], m_valueR[indT][1]);
      vector<double> arrP= linspace(m_valueR[indP][0], m_valueR[indP][2], m_valueR[indP][1]);
      vector<double> arrX= linspace(m_valueR[indX][0], m_valueR[indX][2], m_valueR[indX][1]);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrT.size()*arrP.size()*arrX.size());
      //progressbar
      MultiProgressBar multiBar(arrP.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"3D calculation using "<<m_threadNumOMP<<" threads, T ∈ ["
          <<rangeT[0]<<", "<<rangeT[1]<<"] °C, P ∈ ["
          <<rangeP[0]<<", "<<rangeP[1]<<"] bar, X ∈ ["
          <<rangeX[0]<<", "<<rangeX[1]<<"]"
          <<"\n"<<endl;
      int lenT=arrT.size();
      int lenTX=arrT.size()*arrX.size();
      int lenX = (int)(arrX.size());
      int lenP = (int)(arrP.size());
      #pragma omp parallel for shared(arrT, arrP, arrX, props, lenT, lenTX)
      for (int i = 0; i < lenP; i++)
      {
        for (int j = 0; j < lenX; j++)
        {
          for (int k = 0; k < lenT; k++)
          {
            SWEOS::cH2ONaCl eos;
            eos.prop_pTX(arrP[i]*1e5, arrT[k]+SWEOS::Kelvin, arrX[j]);
            props[k+j*lenT+i*lenTX]=eos.m_prop;
          }
        }
        #pragma omp critical
        multiBar.Update();
      }
      Write2D3DResult(arrT, arrX, arrP, props, m_valueO, "Temperature (°C)", "Salinity", "Pressure (bar)");
    }else if(m_valueV=="PHX" || m_valueV=="PXH" || m_valueV=="HPX" || m_valueV=="HXP" || m_valueV=="XPH" || m_valueV=="XHP")
    {
      int indP=0, indH=1, indX=2;
      if(m_valueV=="PXH")
      {
        indP=0; indX=1; indH=2;
      }else if(m_valueV=="TPX")
      {
        indH=0; indP=1; indX=2;
      }else if(m_valueV=="HXP")
      {
        indH=0; indX=1; indP=2;
      }else if(m_valueV=="XPH")
      {
        indX=0; indP=1; indH=2;
      }else if(m_valueV=="XHP")
      {
        indX=0; indH=1; indP=2;
      }
      double rangeH[2]={m_valueR[indH][0], m_valueR[indH][2]};
      double rangeP[2]={m_valueR[indP][0], m_valueR[indP][2]};
      double rangeX[2]={m_valueR[indX][0], m_valueR[indX][2]};
      double rangePX[4]={m_valueR[indP][0], m_valueR[indP][2], m_valueR[indX][0], m_valueR[indX][2]};
      if(!CheckRanges_H_PX(rangeH[0], rangeH[1], rangePX)) return false;
      if(!CheckRanges_P(rangeP)) return false;
      if(!CheckRanges_X(rangeX)) return false;
      //calculate
      vector<double> arrH= linspace(m_valueR[indH][0], m_valueR[indH][2], m_valueR[indH][1]);
      vector<double> arrP= linspace(m_valueR[indP][0], m_valueR[indP][2], m_valueR[indP][1]);
      vector<double> arrX= linspace(m_valueR[indX][0], m_valueR[indX][2], m_valueR[indX][1]);
      vector<SWEOS::PROP_H2ONaCl> props;
      props.resize(arrH.size()*arrP.size()*arrX.size());
      //progressbar
      MultiProgressBar multiBar(arrP.size(),COLOR_BAR_BLUE);
      omp_set_num_threads(m_threadNumOMP);
      cout<<"3D calculation using "<<m_threadNumOMP<<" threads, H ∈ ["
          <<rangeH[0]<<", "<<rangeH[1]<<"] kJ/kg, P ∈ ["
          <<rangeP[0]<<", "<<rangeP[1]<<"] bar, X ∈ ["
          <<rangeX[0]<<", "<<rangeX[1]<<"]"
          <<"\n"<<endl;
      int lenH=arrH.size();
      int lenHX=arrH.size()*arrX.size();
      int lenX = (int)(arrX.size());
      int lenP = (int)(arrP.size());
      #pragma omp parallel for shared(arrH, arrP, arrX, props, lenH, lenHX)
      for (int i = 0; i < lenP; i++)
      {
        for (int j = 0; j < lenX; j++)
        {
          for (int k = 0; k < lenH; k++)
          {
            SWEOS::cH2ONaCl eos;
            eos.prop_pHX(arrP[i]*1e5, arrH[k]*1000.0, arrX[j]);
            props[k+j*lenH+i*lenHX]=eos.m_prop;
          }
        }
        #pragma omp critical
        multiBar.Update();
      }
      Write2D3DResult(arrH, arrX, arrP, props, m_valueO, "Enthalpy (kJ/kg)", "Salinity", "Pressure (bar)");
    }
    return true;
  }
  bool cSWEOSarg::Validate_1D()
  {
    if(m_valueV.size()!=1)
    {
      cout<<ERROR_COUT<<"if -D set as -D1, the -V option must be one of -VT, -VP, -VX, -VH"<<endl;
      return false;
    }
    if (!m_haveR)
    {
      cout<<ERROR_COUT<<"You set -D1 and -V"<<m_valueV<<", then you must set -R for range of "<<m_valueV<<" in format of -Rmin/delta/max(e.g. -R0/1/100)"<<endl;
      return false;
    }
    if(m_valueR_str.size()!=3)
    {
      cout<<ERROR_COUT<<"You set -D1 and -V"<<m_valueV<<", then you must set -R for range of "<<m_valueV<<" in format of -Rmin/delta/max(e.g. -R0/1/100)"<<endl;
      cout<<ERROR_COUT<<"Option of -R must be 3 values, but what you set is ";
      for(int i=0;i<m_valueR_str.size();i++)cout<<m_valueR_str[i]<<" ";
      cout<<endl;
      return false;
    }
    switch (m_valueV[0])
    {
    case 'T'://change T and fixed P, X
      {
        if(!(m_haveP && CheckRange_P(m_valueP)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument"<<endl;
          return false;
        }
        if(!(m_haveX && CheckRange_X(m_valueX)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument"<<endl;
          return false;
        }
        if(!m_haveO)
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument"<<endl;
          return false;
        }
        double rangeT[2]={m_valueR[0][0], m_valueR[0][2]};
        if(!CheckRanges_T(rangeT)) return false;
        //calculate
        vector<double> arrT= linspace(m_valueR[0][0], m_valueR[0][2], m_valueR[0][1]);
        vector<double> arrP, arrX;
        vector<SWEOS::PROP_H2ONaCl> props;
        SWEOS::cH2ONaCl eos;
        MultiProgressBar multibar(arrT.size(),COLOR_BAR_BLUE);
        for (size_t i = 0; i < arrT.size(); i++)
        {
          arrP.push_back(m_valueP);
          arrX.push_back(m_valueX);
          eos.prop_pTX(arrP[i]*1e5, arrT[i]+SWEOS::Kelvin, arrX[i]);
          props.push_back(eos.m_prop);
          multibar.Update();
        }
        Write1Dresult(m_valueO, arrP, arrX, props);
      }
      break;
    case 'P'://change P and fixed T, X
      {
        if(!(m_haveT && CheckRange_T(m_valueT)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed temperature value by -T argument"<<endl;
          return false;
        }
        if(!(m_haveX && CheckRange_X(m_valueX)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument"<<endl;
          return false;
        }
        if(!m_haveO)
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument"<<endl;
          return false;
        }
        double rangeP[2]={m_valueR[0][0], m_valueR[0][2]};
        if(!CheckRanges_P(rangeP)) return false;
        //calculate
        vector<double> arrP= linspace(m_valueR[0][0], m_valueR[0][2], m_valueR[0][1]);
        vector<double> arrT, arrX;
        vector<SWEOS::PROP_H2ONaCl> props;
        SWEOS::cH2ONaCl eos;
        MultiProgressBar multibar(arrP.size(),COLOR_BAR_BLUE);
        for (size_t i = 0; i < arrP.size(); i++)
        {
          arrT.push_back(m_valueT);
          arrX.push_back(m_valueX);
          eos.prop_pTX(arrP[i]*1e5, arrT[i]+SWEOS::Kelvin, arrX[i]);
          props.push_back(eos.m_prop);
          multibar.Update();
        }
        Write1Dresult(m_valueO, arrP, arrX, props);
      }
      break;
    case 'X'://change X and fixed P, T
      {
        if(!(m_haveT && CheckRange_T(m_valueP)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed temperature value by -T argument"<<endl;
          return false;
        }
        if(!(m_haveP && CheckRange_P(m_valueP)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument"<<endl;
          return false;
        }
        if(!m_haveO)
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument"<<endl;
          return false;
        }
        double rangeX[2]={m_valueR[0][0], m_valueR[0][2]};
        if(!CheckRanges_X(rangeX)) return false;
        //calculate
        vector<double> arrX= linspace(m_valueR[0][0], m_valueR[0][2], m_valueR[0][1]);
        vector<double> arrT, arrP;
        vector<SWEOS::PROP_H2ONaCl> props;
        SWEOS::cH2ONaCl eos;
        MultiProgressBar multibar(arrX.size(),COLOR_BAR_BLUE);
        for (size_t i = 0; i < arrX.size(); i++)
        {
          arrT.push_back(m_valueT);
          arrP.push_back(m_valueP);
          eos.prop_pTX(arrP[i]*1e5, arrT[i]+SWEOS::Kelvin, arrX[i]);
          props.push_back(eos.m_prop);
          multibar.Update();
        }
        Write1Dresult(m_valueO, arrP, arrX, props);
      }
      break;
    case 'H'://change H and fixed P, X
      {
        if(!(m_haveP && CheckRange_P(m_valueP)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed pressure value by -P argument"<<endl;
          return false;
        }
        if(!(m_haveX && CheckRange_X(m_valueX)))
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", but you don't set a proper fixed salinity value by -X argument"<<endl;
          return false;
        }
        if(!m_haveO)
        {
          cout<<ERROR_COUT<<"Selected calculation mode is 1D calculation, change "<<m_valueV<<", you have to specify an output file by -O argument"<<endl;
          return false;
        }
        // if(!checkranges)
        double rangeH[2]={m_valueR[0][0], m_valueR[0][2]};
        if(!CheckRanges_H(rangeH, m_valueP, m_valueX)) return false;
        //calculate
        vector<double> arrH= linspace(m_valueR[0][0], m_valueR[0][2], m_valueR[0][1]);
        vector<double> arrP, arrX;
        vector<SWEOS::PROP_H2ONaCl> props;
        SWEOS::cH2ONaCl eos;
        MultiProgressBar multibar(arrH.size(),COLOR_BAR_BLUE);
        for (size_t i = 0; i < arrH.size(); i++)
        {
          arrP.push_back(m_valueP);
          arrX.push_back(m_valueX);
          eos.prop_pHX(arrP[i]*1e5, arrH[i]*1000.0, arrX[i]);
          props.push_back(eos.m_prop);
          multibar.Update();
        }
        Write1Dresult(m_valueO, arrP, arrX, props);
      }
      break;
    default:
      break;
    }
    
    return true;
  }
  template<typename T>
  vector<T> cSWEOSarg::linspace(T xmin, T xmax, T dx)
  {
    vector<T> tmp;
    for (T x = xmin; x <= xmax; x=x+dx)
    {
      tmp.push_back(x);
    }
    return tmp;
  }
  bool cSWEOSarg::Validate_0D()
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
        if(!CheckRange_P(m_valueP))return false;
        if(!CheckRange_T(m_valueT))return false;
        if(!CheckRange_X(m_valueX))return false;
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
          cout<<WARN_COUT<<"selected calculation mode is single point PHX, you must specify enthalpy by -H"<<endl;
        }
        if(!m_haveP)
        {
          cout<<WARN_COUT<<"selected calculation mode is single point PHX, you must specify pressure by -P"<<endl;
        }
        if(!m_haveX)
        {
          cout<<WARN_COUT<<"selected calculation mode is single point PHX, you must specify salinity by -X"<<endl;
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
          <<COLOR_PURPLE<<"\tPlease make sure your input file with 3 columns in order of "<<COLOR_RED<<m_valueV<<COLOR_DEFAULT<<endl;
          calculateMultiPoints_PTX_PHX(m_valueV,m_valueG, m_valueO,"H");
        }
        else
        {
          cout<<ERROR_COUT<<"There neither full -H, -P, -X options nor -G argument, swEOS will exit"<<endl;
          exit(0);
        } 
    }
    else
    {
      cout<<ERROR_COUT<<"unknown -D parameter: "<<m_valueV<<endl;
      return false;
    }
      return true;
  }
  bool cSWEOSarg::CheckRange_T(double T0, double TMIN, double TMAX)
  {
    if(T0>SWEOS::TMAX-SWEOS::Kelvin || T0<SWEOS::TMIN-SWEOS::Kelvin)
    {
      cout<<ERROR_COUT<<"-T specify temperature ="<<T0<<"°C = "<<T0+SWEOS::Kelvin<<" K, out of range ["<<TMIN<<", "<<TMAX<<"] K"<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_T(double Trange[2], double TMIN, double TMAX)
  {
    if(Trange[0]>SWEOS::TMAX-SWEOS::Kelvin || Trange[0]<SWEOS::TMIN-SWEOS::Kelvin)
    {
      cout<<ERROR_COUT<<"The minimum value of temperature specified by -R argument Tmin ="<<Trange[0]<<"°C = "<<Trange[0]+SWEOS::Kelvin<<" K, out of range ["<<TMIN<<", "<<TMAX<<"] K"<<endl;
      return false;
    }
    if(Trange[1]>SWEOS::TMAX-SWEOS::Kelvin || Trange[1]<SWEOS::TMIN-SWEOS::Kelvin)
    {
      cout<<ERROR_COUT<<"The maximum value of temperature specified by -R argument Tmax ="<<Trange[1]<<"°C = "<<Trange[1]+SWEOS::Kelvin<<" K, out of range ["<<TMIN<<", "<<TMAX<<"] K"<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRange_P(double P0, double PMIN, double PMAX)
  {
    if(P0>PMAX/1e5 || P0<PMIN/1e5)
    {
      cout<<ERROR_COUT<<"-P specify pressure ="<<P0<<"bar = "<<P0*1e5<<" Pa, out of range ["<<PMIN<<", "<<PMAX<<"] Pa"<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_P(double Prange[2], double PMIN, double PMAX)
  {
    if(Prange[0]>PMAX/1e5 || Prange[0]<PMIN/1e5)
    {
      cout<<ERROR_COUT<<"The minimum value of pressure specified by -R argument Pmin ="<<Prange[0]<<"bar = "<<Prange[0]*1e5<<" Pa, out of range ["<<PMIN<<", "<<PMAX<<"] Pa"<<endl;
      return false;
    }
    if(Prange[1]>PMAX/1e5 || Prange[1]<PMIN/1e5)
    {
      cout<<ERROR_COUT<<"The minimum value of pressure specified by -R argument Pmax ="<<Prange[1]<<"bar = "<<Prange[1]*1e5<<" Pa, out of range ["<<PMIN<<", "<<PMAX<<"] Pa"<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_X(double Xrange[2], double XMIN, double XMAX)
  {
    if(Xrange[0]>XMAX || Xrange[0]<XMIN)
    {
      cout<<ERROR_COUT<<"The minimum value of salinity specified by -R argument Xmin ="<<Xrange[0]<<", out of range ["<<XMIN<<", "<<XMAX<<"]"<<endl;
      return false;
    }
    if(Xrange[1]>XMAX || Xrange[1]<XMIN)
    {
      cout<<ERROR_COUT<<"The maximum value of salinity specified by -R argument Xmax ="<<Xrange[1]<<", out of range ["<<XMIN<<", "<<XMAX<<"]"<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRange_X(double X0, double XMIN, double XMAX)
  {
    if(X0>XMAX || X0<XMIN)
    {
      cout<<ERROR_COUT<<"-X specify salinity ="<<X0<<" wt. % NaCl, out of range ["<<XMIN<<", "<<XMAX<<"] wt. % NaCl"<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRange_H(double H0, double P0, double X0)
  {
    double HMIN=1e30, HMAX=-1e30;
    double Trange[2]={SWEOS::TMIN, SWEOS::TMAX};
    double P,T,X;
    SWEOS::cH2ONaCl eos;
    // for (int i = 0; i < 2; i++)
    {
      P=P0;
      for (int j = 0; j < 2; j++)
      {
        T=Trange[j];
        // for (int k = 0; k < 2; k++)
        {
          X=X0;
          eos.prop_pTX(P*1e5, T, X);
          HMIN=(eos.m_prop.H<HMIN ? eos.m_prop.H : HMIN);
          HMAX=(eos.m_prop.H>HMAX ? eos.m_prop.H : HMAX);
        }
      }
    }
    if(H0<HMIN/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The enthalpy value is specified by -H argument is H="<<H0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the minimum enthalpy corresponding to P="<<P0<<" bar, X="<<X0<<" is "<<HMIN/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    if(H0>HMAX/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The enthalpy value is specified by -H argument is H="<<H0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the maximum enthalpy corresponding to P="<<P0<<" bar, X="<<X0<<" is "<<HMAX/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_H(double Hrange[2], double P0, double X0)
  {
    double HMIN=1e30, HMAX=-1e30;
    double Trange[2]={SWEOS::TMIN, SWEOS::TMAX};
    double P,T,X;
    SWEOS::cH2ONaCl eos;
    // for (int i = 0; i < 2; i++)
    {
      P=P0;
      for (int j = 0; j < 2; j++)
      {
        T=Trange[j];
        // for (int k = 0; k < 2; k++)
        {
          X=X0;
          eos.prop_pTX(P*1e5, T, X);
          HMIN=(eos.m_prop.H<HMIN ? eos.m_prop.H : HMIN);
          HMAX=(eos.m_prop.H>HMAX ? eos.m_prop.H : HMAX);
        }
      }
    }

    if(Hrange[0]<HMIN/1000 || Hrange[0]>HMAX/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The minimum value of enthalpy is specified by -R argument is Hmin="<<Hrange[0]<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the range of enthalpy corresponding to P="<<P0<<" bar, X="<<X0<<" is ["<<HMIN/1000<<", "<<HMAX/1000<<"] kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    if(Hrange[1]<HMIN/1000 || Hrange[1]>HMAX/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The maximum value of enthalpy is specified by -R argument is Hmin="<<Hrange[1]<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the range of enthalpy corresponding to P="<<P0<<" bar, X="<<X0<<" is ["<<HMIN/1000<<", "<<HMAX/1000<<"] kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_H_PX(double HMIN0, double HMAX0, double PXrange[4])
  {
    double HMIN=1e30, HMAX=-1e30;
    double Trange[2]={SWEOS::TMIN, SWEOS::TMAX};
    double P,T,X;
    SWEOS::cH2ONaCl eos;
    for (int i = 0; i < 2; i++)
    {
      P=PXrange[i];
      for (int j = 0; j < 2; j++)
      {
        T=Trange[j];
        for (int k = 0; k < 2; k++)
        {
          X=PXrange[k+2];
          eos.prop_pTX(P*1e5, T, X);
          HMIN=(eos.m_prop.H<HMIN ? eos.m_prop.H : HMIN);
          HMAX=(eos.m_prop.H>HMAX ? eos.m_prop.H : HMAX);
        }
      }
    }
    if(HMIN0<HMIN/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The minimum value of enthalpy is specified by -R argument is MIN="<<HMIN0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the minimum enthalpy corresponding to P["<<PXrange[0]<<", "<<PXrange[1]<<"] bar, X["<<PXrange[2]<<", "<<PXrange[3]<<"] is "<<HMIN/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    if(HMAX0>HMAX/1000) 
    {
      cout<<WARN_COUT<<COLOR_RED<<"The maximum value of enthalpy is specified by -R argument is MAX="<<HMAX0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the maximum enthalpy corresponding to P["<<PXrange[0]<<", "<<PXrange[1]<<"] bar, X["<<PXrange[2]<<", "<<PXrange[3]<<"] is "<<HMAX/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_H_P(double HMIN0, double HMAX0, double Prange[2], double X0)
  {
    double HMIN=1e30, HMAX=-1e30;
    double Trange[2]={SWEOS::TMIN, SWEOS::TMAX};
    double P,T,X;
    SWEOS::cH2ONaCl eos;
    for (int i = 0; i < 2; i++)
    {
      P=Prange[i];
      for (int j = 0; j < 2; j++)
      {
        T=Trange[j];
        // for (int k = 0; k < 2; k++)
        {
          X=X0;
          eos.prop_pTX(P*1e5, T, X);
          HMIN=(eos.m_prop.H<HMIN ? eos.m_prop.H : HMIN);
          HMAX=(eos.m_prop.H>HMAX ? eos.m_prop.H : HMAX);
        }
      }
    }
    if(HMIN0<HMIN/1000)
    {
      cout<<ERROR_COUT<<COLOR_RED<<"The minimum value of enthalpy is specified by -R argument is MIN="<<HMIN0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the minimum enthalpy corresponding to P["<<Prange[0]<<", "<<Prange[1]<<"] bar, X="<<X0<<" is "<<HMIN/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    if(HMAX0>HMAX/1000)
    {
      cout<<ERROR_COUT<<COLOR_RED<<"The maximum value of enthalpy is specified by -R argument is MAX="<<HMAX0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the maximum enthalpy corresponding to P["<<Prange[0]<<", "<<Prange[1]<<"] bar, X="<<X0<<"] is "<<HMAX/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    return true;
  }
  bool cSWEOSarg::CheckRanges_H_X(double HMIN0, double HMAX0, double Xrange[2], double P0)
  {
    double HMIN=1e30, HMAX=-1e30;
    double Trange[2]={SWEOS::TMIN, SWEOS::TMAX};
    double P,T,X;
    SWEOS::cH2ONaCl eos;
    // for (int i = 0; i < 2; i++)
    {
      P=P0;
      for (int j = 0; j < 2; j++)
      {
        T=Trange[j];
        for (int k = 0; k < 2; k++)
        {
          X=Xrange[k];
          eos.prop_pTX(P*1e5, T, X);
          HMIN=(eos.m_prop.H<HMIN ? eos.m_prop.H : HMIN);
          HMAX=(eos.m_prop.H>HMAX ? eos.m_prop.H : HMAX);
        }
      }
    }
    if(HMIN0<HMIN/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The minimum value of enthalpy is specified by -R argument is MIN="<<HMIN0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the minimum enthalpy corresponding to P="<<P0<<" bar, X["<<Xrange[0]<<", "<<Xrange[1]<<"] is "<<HMIN/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    if(HMAX0>HMAX/1000)
    {
      cout<<WARN_COUT<<COLOR_RED<<"The maximum value of enthalpy is specified by -R argument is MAX="<<HMAX0<<" kJ/kg"
          <<"may be out of range and could cause swEOS crash.\n"
          <<"Because the maximum enthalpy corresponding to P="<<P0<<" bar, X["<<Xrange[0]<<", "<<Xrange[1]<<"] is "<<HMAX/1000<<" kJ/kg"
          <<COLOR_DEFAULT<<endl;
      return false;
    }
    return true;
  }
  vector<string> string_split (string s, string delimiter) {
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
  vector<SWEOS::PROP_H2ONaCl> calculateMultiPoints_PTX_PHX(string valueV, string filePTX, string outFile, string isT_H)
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

    //calculate
    SWEOS::cH2ONaCl eos;
    vector<SWEOS::PROP_H2ONaCl> props;
    MultiProgressBar multibar(P.size(),COLOR_BAR_BLUE);
    if(isT_H=="T")
    {
      for(int i=0;i<P.size();i++)
      {
        eos.prop_pTX(P[i]*1e5, T_H[i]+SWEOS::Kelvin, X[i]);
        props.push_back(eos.m_prop);
        multibar.Update();
      }
    }else if(isT_H=="H")
    {
      for(int i=0;i<P.size();i++)
      {
        eos.prop_pHX(P[i]*1e5, T_H[i]*1000, X[i]);
        props.push_back(eos.m_prop);
        multibar.Update();
      }
    }
    //write to file
    Write1Dresult(outFile,P,X,props);

    return props;
  }
  bool Write1Dresult(string outFile,vector<double> P, vector<double> X, vector<SWEOS::PROP_H2ONaCl> props)
  {
    string extname;
    vector<string> tmp=string_split(outFile,".");
    if(tmp.size()>=1)
    {
      extname=tmp[tmp.size()-1];
    }
    if(extname=="csv")
    {
      WriteCSV(outFile, P, X, props);
    }else
    {
      cout<<WARN_COUT<<"Unrecognized format: "<<outFile<<endl;
      cout<<COLOR_GREEN<<"Write results into csv file format"<<COLOR_DEFAULT<<endl;
      WriteCSV(outFile, P, X, props);
    }
    return true;
  }
  bool WriteCSV(string outFile,vector<double> P, vector<double> X, vector<SWEOS::PROP_H2ONaCl> props)
  {
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
      for(int i=0;i<P.size();i++)
      {
        fout<<props[i].T<<", "<<P[i]<<", "<<X[i]<<", "<<props[i].Region<<", "<<eos.m_phaseRegion_name[props[i].Region]<<", "
            <<props[i].Rho<<", "<<props[i].Rho_l<<", "<<props[i].Rho_v<<", "<<props[i].Rho_h<<", "
            <<props[i].H/1000.0<<", "<<props[i].H_l/1000.0<<", "<<props[i].H_v/1000.0<<", "<<props[i].H_h/1000.0<<", "
            <<props[i].S_l<<", "<<props[i].S_v<<", "<<props[i].S_h<<", "
            <<props[i].Mu_l<<", "<<props[i].Mu_v<<", "
            <<props[i].X_l<<", "<<props[i].X_v
            <<endl;
      }
      fout.close();
      cout<<COLOR_BLUE<<"Results have been saved to file: "<<outFile<<endl;
    }else
    {
      if(outFile=="")
      {
        cout<<WARN_COUT<<"The output file has to be specified by -O argument"<<outFile<<endl;
        exit(0);
      }

      if(!fout)cout<<WARN_COUT<<"The output file name is specified by -O argument, but open file failed: "<<outFile
                   <<endl;
      exit(0);
      // cout<<"T(C), P(bar), X, Phase Index, Phase Region, "
      //     <<"Bulk density(kg/m3), Liquid density(kg/m3), Vapour density(kg/m3), Halite density(kg/m3), "
      //     <<"Bulk enthalpy(kJ/kg), Liquid enthalpy(kJ/kg), Vapour enthalpy(kJ/kg), Halite enthalpy(kJ/kg), "
      //     <<"Liquid Saturation, Vapour Saturation, Halite Saturation, "
      //     <<"Liquid viscosity, Vapour viscosity, "
      //     <<"Liquid salinity, Vapour salinity"
      //     <<endl;
      // for(int i=0;i<P.size();i++)
      // {
      //   cout<<props[i].T+SWEOS::Kelvin<<", "<<P[i]<<", "<<X[i]<<", "<<props[i].Region<<", "<<eos.m_phaseRegion_name[props[i].Region]<<", "
      //       <<props[i].Rho<<", "<<props[i].Rho_l<<", "<<props[i].Rho_v<<", "<<props[i].Rho_h<<", "
      //       <<props[i].H<<", "<<props[i].H_l<<", "<<props[i].H_v<<", "<<props[i].H_h<<", "
      //       <<props[i].S_l<<", "<<props[i].S_v<<", "<<props[i].S_h<<", "
      //       <<props[i].Mu_l<<", "<<props[i].Mu_v<<", "
      //       <<props[i].X_l<<", "<<props[i].X_v
      //       <<endl;
      // }
    }
    return true;
  }
}