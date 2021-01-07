#include<iostream>
#include "Polynomial.h"
using namespace std;

#include"H2ONaCl.H"
using namespace H2ONaCl;

// void PhaseRegion_PTX(double P, double T, double X);
// void Prop_PHX(double P, double H_MJ, double Xwt);
// void PhaseRegion3D_PTX();
// void Figure2_Driesner2007(double T, string filename, double pmax=400e5);
// void Figure2_Driesner2007_log(double T, string filename, double pmax=400e5, double xmin=-11);
// void test_root_poly(double P)
// {
//     const int N = 10;
//     double f[11] = {0.00464, 5E-07, 16.9078, -269.148, 7632.04, -49563.6, 233119.0, -513556.0, 549708.0, -284628.0, NaCl::P_Triple};
//     double P_VLH = 0;
//     for (size_t i = 0; i < 10; i++) //calculate P_VLH (eq. 10) and f[10] simultaneously
//     {
//         f[10] -= f[i];
//     }

//     Polynomial polynomial;
//     int degree = 10;
//     std::vector<double> coefficient_vector;
//     coefficient_vector.resize(degree + 1);
//     double * coefficient_vector_ptr = &coefficient_vector[0];
//     for (size_t i = 0; i < degree+1; i++)
//     {
//         coefficient_vector_ptr[i]=f[i];
//     }
//     coefficient_vector_ptr[0]=coefficient_vector_ptr[0]-P;
//     polynomial.SetCoefficients(coefficient_vector_ptr, degree);

//     std::vector<double> real_vector;
//     std::vector<double> imag_vector;

//     // int degree = polynomial.Degree();

//     real_vector.resize(degree);
//     imag_vector.resize(degree);

//     double * real_vector_ptr = &real_vector[0];
//     double * imag_vector_ptr = &imag_vector[0];

//     int root_count= 0;

//     if (polynomial.FindRoots(real_vector_ptr,
//                                 imag_vector_ptr,
//                                 &root_count) == PolynomialRootFinder::SUCCESS)
//     {
//         int i = 0;

//         for (i = 0; i < root_count; ++i)
//         {
//             if(imag_vector_ptr[i]==0)
//             std::cout << "Root " << i << " = " << real_vector_ptr[i]*NaCl::T_Triple << " + i " << imag_vector_ptr[i] << std::endl;
//         }
//     }
//     else
//     {
//         std::cout << "Failed to find all roots." << std::endl;
//     }
// }
int main()
{
    // test_root_poly(350);
    // 1. one point calculation
    // PhaseRegion_PTX(31600000,100+273.15,0.3);
    // Prop_PHX(500e5, 3e6, 0.03);

    // 2. 3D P-T-X calculation
    // PhaseRegion3D_PTX();

    // 3. Benchmark test and compar with Driesner 2007
    // Figure2_Driesner2007(300, "Figure2a_Driesner2007.csv",400e5);
    // Figure2_Driesner2007(373.976, "Figure2c_Driesner2007.csv",400e5);
    // Figure2_Driesner2007(375, "Figure2e_Driesner2007.csv",400e5);
    // Figure2_Driesner2007(500, "Figure2g_Driesner2007.csv",800e5);
    // Figure2_Driesner2007(800, "Figure2i_Driesner2007.csv",2500e5);
    // Figure2_Driesner2007(1000, "Figure2k_Driesner2007.csv",2500e5);

    // Figure2_Driesner2007_log(300, "Figure2b_Driesner2007.csv",400e5);
    // Figure2_Driesner2007_log(373.976, "Figure2d_Driesner2007.csv",400e5);
    // Figure2_Driesner2007_log(375, "Figure2f_Driesner2007.csv",400e5);
    // Figure2_Driesner2007_log(500, "Figure2h_Driesner2007.csv",800e5);
    // Figure2_Driesner2007_log(800, "Figure2j_Driesner2007.csv",100e5,-3);

    return 0;
}
// void Prop_PHX(double P, double H, double Xwt)
// {
//     cout.precision(8);//control cout precision of float
//     cH2ONaCl eos;
//     // H2ONaCl::PROP_H2ONaCl prop;
//     eos.m_prop=eos.prop_pHX(P,H,Xwt);
//     eos.setColorPrint(true);
//     cout<<eos<<endl;
//     cout<<"Region: "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
//     cout<<"T: "<<eos.m_prop.T<<" Rho: "<<eos.m_prop.Rho<<" H: "<<eos.m_prop.H<<endl;
//     cout<<"Rho_l: "<<eos.m_prop.Rho_l<<" Rho_v: "<<eos.m_prop.Rho_v<<" Rho_h: "<<eos.m_prop.Rho_h<<endl;
//     cout<<"h_l: "<<eos.m_prop.H_l<<" h_v: "<<eos.m_prop.H_v<<" h_h: "<<eos.m_prop.H_h<<endl;
//     cout<<"Rho: "<<eos.m_prop.Rho<<" H: "<<eos.m_prop.H<<endl;
//     cout<<"S_l: "<<eos.m_prop.S_l<<" S_v: "<<eos.m_prop.S_v<<" S_h: "<<eos.m_prop.S_h<<endl;
//     cout<<"Mu_l: "<<eos.m_prop.Mu_l<<" Mu_v: "<<eos.m_prop.Mu_v<<endl;
// }

// void Figure2_Driesner2007_log(double T, string filename,double pmax, double xmin)
// {
//     double xmax=2, dx=0.01;
//     double pmin=1e5,dp=1e5;
//     ofstream fpout(filename);
//     if(!fpout)
//     {
//         cout<<"Error: File can not opened->Figure2_Driesner2007a: "<<filename<<endl;
//         exit(0);
//     }
//     fpout<<"T,P,X,Region,Rho,H,Rho_l,Rho_v,Rho_h,Mu_l,Mu_v"<<endl;
//     for (double logX = xmin; logX < xmax; logX=logX+dx)
//     {
//         double X=pow(10,logX)/100.0;
//         for (double P = pmin; P < pmax; P=P+dp)
//         {
//             cH2ONaCl eos;
//             eos.m_prop=eos.prop_pTX(P,T,X);
//             fpout<<T<<", "<<P<<", "<<logX<<", "<<eos.m_prop.Region<<", "
//                  <<eos.m_prop.Rho<<", "<<eos.m_prop.H<<", "
//                  <<eos.m_prop.Rho_l<<", "<<eos.m_prop.Rho_v<<", "<<eos.m_prop.Rho_h<<", "
//                  <<eos.m_prop.Mu_l<<", "<<eos.m_prop.Mu_v<<endl;
//         }
//     }
//     fpout.close();
// }
// void Figure2_Driesner2007(double T, string filename,double pmax)
// {
//     double xmin=0, xmax=1, dx=0.001;
//     double pmin=1e5,dp=1e5;
//     ofstream fpout(filename);
//     if(!fpout)
//     {
//         cout<<"Error: File can not opened->Figure2_Driesner2007a: "<<filename<<endl;
//         exit(0);
//     }
//     fpout<<"T,P,X,Region,Rho,H,Rho_l,Rho_v,Rho_h,Mu_l,Mu_v"<<endl;
//     for (double X = xmin; X < xmax; X=X+dx)
//     {
//         for (double P = pmin; P < pmax; P=P+dp)
//         {
//             cH2ONaCl eos;
//             eos.m_prop=eos.prop_pTX(P,T,X);
//             fpout<<T<<", "<<P<<", "<<X<<", "<<eos.m_prop.Region<<", "
//                  <<eos.m_prop.Rho<<", "<<eos.m_prop.H<<", "
//                  <<eos.m_prop.Rho_l<<", "<<eos.m_prop.Rho_v<<", "<<eos.m_prop.Rho_h<<", "
//                  <<eos.m_prop.Mu_l<<", "<<eos.m_prop.Mu_v<<endl;
//         }
//     }
//     fpout.close();
// }

// void PhaseRegion_PTX(double P, double T, double X)
// {
//     cout.precision(8);//control cout precision of float
//     cH2ONaCl eos;
//     eos.m_prop=eos.prop_pTX(P,T,X);

//     cout<<"Region: "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
//     cout<<"Rho: "<<eos.m_prop.Rho<<" H: "<<eos.m_prop.H<<endl;
//     cout<<"Rho_l: "<<eos.m_prop.Rho_l<<" Rho_v: "<<eos.m_prop.Rho_v<<" Rho_h: "<<eos.m_prop.Rho_h<<endl;
//     cout<<"h_l: "<<eos.m_prop.H_l<<" h_v: "<<eos.m_prop.H_v<<" h_h: "<<eos.m_prop.H_h<<endl;
//     cout<<"Rho: "<<eos.m_prop.Rho<<" H: "<<eos.m_prop.H<<endl;
//     cout<<"S_l: "<<eos.m_prop.S_l<<" S_v: "<<eos.m_prop.S_v<<" S_h: "<<eos.m_prop.S_h<<endl;
//     cout<<"Mu_l: "<<eos.m_prop.Mu_l<<" Mu_v: "<<eos.m_prop.Mu_v<<endl;
// }
// void PhaseRegion3D_PTX()
// {
//     double dT=10;
//     double dP=10e5;
//     double dX=0.02;
//     double MAXT=1000;
//     double MINT=dT;
//     double MAXP=2200E5;
//     double MINP=1E5;
//     double MINX=0;
//     double MAXX=1;

//     vector<double> T, P, X;
//     for(double t=MINT;t<=MAXT;t=t+dT)
//     {
//         T.push_back(t);
//     }
//     for(double p=MINP;p<=MAXP;p=p+dP)
//     {
//         P.push_back(p);
//     }
//     for(double x=MINX;x<=MAXX;x=x+dX)
//     {
//         X.push_back(x);
//     }
//     vector<PROP_H2ONaCl> props;

//     // ofstream fpout("region.csv");
//     // fpout.precision(8);//control cout precision of float
//     // if(!fpout)
//     // {
//     //     cout<<"ERROR: Can not open file: region.csv"<<endl;
//     //     exit(0);
//     // }
//     // for (size_t i = 0; i < X.size(); i++)
//     // {
//     //     for (size_t j = 0; j < T.size(); j++)
//     //     {
//     //             for (size_t k = 0; k < P.size(); k++)
//     //             {
//     //                 cH2ONaCl eos(P[k],T[j],X[i]);
//     //                 eos.Calculate(); 
//     //                 fpout<<T[j]<<","<<P[k]<<","<<X[i]<<","<<eos.m_phaseRegion_name[eos.m_prop.Region]
//     //                      <<","<<eos.m_prop.X_l<<","<<eos.m_prop.X_v<<endl;
//     //             }
//     //     }
//     //     cout<<i<<endl;
//     // }
//     // fpout.close();

//     // write to VTK
//     for (size_t i = 0; i < X.size(); i++)
//     {
//         for (size_t j = 0; j < P.size(); j++)
//         {
//                 for (size_t k = 0; k < T.size(); k++)
//                 {
//                     cH2ONaCl eos;
//                     eos.m_prop=eos.prop_pTX(P[j],T[k],X[i]);
//                     props.push_back(eos.m_prop);
//                     // cout<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
//                     // fpout<<eos.m_prop.Region<<" ";
//                 }
//                 // fpout<<endl;
//         }
//         cout<<i<<endl;
//     }
    
//     // fpout.close();
//     cH2ONaCl eos;
//     eos.writeProps2VTK(T,P,X,props,"region.vtk");
// }