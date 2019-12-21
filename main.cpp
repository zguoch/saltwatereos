#include<iostream>

using namespace std;

#include"H2ONaCl.h"
void PhaseRegion_PTX(double P, double T, double X);
void PhaseRegion3D_PTX();
int main()
{
    // 1. one point calculation
    // PhaseRegion_PTX(43000000,155,0.3);

    // 2. 3D P-T-X calculation
    PhaseRegion3D_PTX();
    return 0;
}
void PhaseRegion_PTX(double P, double T, double X)
{
    cout.precision(8);//control cout precision of float
    cH2ONaCl eos(P,T,X);
    eos.Calculate();

    cout<<"Region: "<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
}
void PhaseRegion3D_PTX()
{
    double dT=1;
    double dP=1e5;
    double dX=0.1;
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

    // ofstream fpout("region.csv");
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
    //                 fpout<<T[j]<<","<<P[k]<<","<<X[i]<<","<<eos.m_phaseRegion_name[eos.m_prop.Region]<<endl;
    //             }
    //     }
    //     cout<<i<<endl;
    // }

    // fpout.close();


    // calculate and write to file
    ofstream fpout("region.vtk");
    if(!fpout)
    {
        cout<<"ERROR: Can not open file: region.vtk"<<endl;
        exit(0);
    }
    //  write vtk head
    fpout<<"# vtk DataFile Version 2.0"<<endl;
    fpout<<"Properties of seawater"<<endl;
    fpout<<"ASCII"<<endl;
    fpout<<"DATASET RECTILINEAR_GRID"<<endl;
    fpout<<"DIMENSIONS "<<T.size()<<" "<<P.size()<<" "<<X.size()<<endl;
    fpout<<"X_COORDINATES "<<T.size()<<" float"<<endl;
    for(int i=0;i<T.size();i++)fpout<<(T[i]-TMIN)/(TMAX-TMIN)<<" ";fpout<<endl;
    fpout<<"Y_COORDINATES "<<P.size()<<" float"<<endl;
    for(int i=0;i<P.size();i++)fpout<<(P[i]-PMIN)/(PMAX-PMIN)<<" ";fpout<<endl;
    fpout<<"Z_COORDINATES "<<X.size()<<" float"<<endl;
    for(int i=0;i<X.size();i++)fpout<<(X[i]-XMIN)/(XMAX-XMIN)<<" ";fpout<<endl;

    fpout<<"POINT_DATA "<<T.size()*P.size()*X.size()<<endl;
    fpout<<"SCALARS region float"<<endl;
    fpout<<"LOOKUP_TABLE default"<<endl;
    for (size_t i = 0; i < X.size(); i++)
    {
        for (size_t j = 0; j < P.size(); j++)
        {
                for (size_t k = 0; k < T.size(); k++)
                {
                    cH2ONaCl eos(P[j],T[k],X[i]);
                    eos.Calculate(); 
                    // cout<<eos.m_phaseRegion_name[eos.m_prop.Region]<<" Xl: "<<eos.m_prop.X_l<<" Xv: "<<eos.m_prop.X_v<<endl;
                    fpout<<eos.m_prop.Region<<" ";
                }
                fpout<<endl;
        }
        cout<<i<<endl;
    }
    
    fpout.close();
}