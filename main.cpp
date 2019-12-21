#include<iostream>

using namespace std;

#include"H2ONaCl.h"
void PhaseRegion_PTX(double P, double T, double X);
void PhaseRegion3D_PTX();
int main()
{
    // 1. one point calculation
    PhaseRegion_PTX(31600000,100,0.3);

    // 2. 3D P-T-X calculation
    // PhaseRegion3D_PTX();
    return 0;
}
void PhaseRegion_PTX(double P, double T, double X)
{
    cout.precision(8);//control cout precision of float
    cH2ONaCl eos(P,T,X);
    eos.Calculate();

    cout<<"Region: "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
    cout<<"Rho_l: "<<eos.m_prop.Rho_l<<" Rho_v: "<<eos.m_prop.Rho_v<<" Rho_h: "<<eos.m_prop.Rho_h<<endl;
    cout<<"h_l: "<<eos.m_prop.H_l<<" h_v: "<<eos.m_prop.H_v<<" h_h: "<<eos.m_prop.H_h<<endl;
    cout<<"Rho: "<<eos.m_prop.Rho<<" H: "<<eos.m_prop.H<<endl;
    cout<<"S_l: "<<eos.m_prop.S_l<<" S_v: "<<eos.m_prop.S_v<<" S_h: "<<eos.m_prop.S_h<<endl;
    cout<<"Mu_l: "<<eos.m_prop.Mu_l<<" Mu_v: "<<eos.m_prop.Mu_v<<endl;
}
void PhaseRegion3D_PTX()
{
    double dT=10;
    double dP=10e5;
    double dX=0.02;
    double MAXT=500;
    double MINT=dT;
    double MAXP=500E5;
    double MINP=1E5;
    double MINX=0;
    double MAXX=1;

    vector<double> T, P, X;
    for(double t=MINT;t<=MAXT;t=t+dT)
    {
        T.push_back(t);
    }
    for(double p=MINP;p<=MAXP;p=p+dP)
    {
        P.push_back(p);
    }
    for(double x=MINX;x<=MAXX;x=x+dX)
    {
        X.push_back(x);
    }
    vector<PROP_H2ONaCl> props;

    // ofstream fpout("region.csv");
    // fpout.precision(8);//control cout precision of float
    // if(!fpout)
    // {
    //     cout<<"ERROR: Can not open file: region.csv"<<endl;
    //     exit(0);
    // }
    // for (size_t i = 0; i < X.size(); i++)
    // {
    //     for (size_t j = 0; j < T.size(); j++)
    //     {
    //             for (size_t k = 0; k < P.size(); k++)
    //             {
    //                 cH2ONaCl eos(P[k],T[j],X[i]);
    //                 eos.Calculate(); 
    //                 fpout<<T[j]<<","<<P[k]<<","<<X[i]<<","<<eos.m_phaseRegion_name[eos.m_prop.Region]
    //                      <<","<<eos.m_prop.X_l<<","<<eos.m_prop.X_v<<endl;
    //             }
    //     }
    //     cout<<i<<endl;
    // }
    // fpout.close();

    // write to VTK
    for (size_t i = 0; i < X.size(); i++)
    {
        for (size_t j = 0; j < P.size(); j++)
        {
                for (size_t k = 0; k < T.size(); k++)
                {
                    cH2ONaCl eos(P[j],T[k],X[i]);
                    eos.Calculate();
                    props.push_back(eos.m_prop);
                    // cout<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
                    // fpout<<eos.m_prop.Region<<" ";
                }
                // fpout<<endl;
        }
        cout<<i<<endl;
    }
    
    // fpout.close();
    cH2ONaCl eos;
    eos.writeProps2VTK(T,P,X,props,"region.vtk");
}