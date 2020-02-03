/**
 * @file IAPWS-IF97.H
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief 
 * @version 1.0
 * @date 2019-09-27
 * 
 * Using freesteam, you can set the state of steam using various pairs of know properties. The following input pairs are implemented:
	(p,T)
	(p,h)
	(p,s)
	(p,v)
	(T,s)
	(T,x)
In C, these states are set using syntax like SteamState S = freesteam_set_ph(p,h). In Python, the syntax is S = steam_ph(p,h). The routines return a 'SteamState' object that can then be queried for any output properties you wish to know, including
	p	Pressure
	T	Temperature
	h	Enthalpy
	s	Entropy
	cp	Isobaric heat capacity
	cv	Isochoric heat capacity
	v	Specific volume
	rho	Density
	u	Internal energy
	x Quality (when saturated)
	k	Thermal conductivity
	mu	Dynamic viscosity
	w	Speed of sound

 * @copyright Copyright (c) 2019
 * 
 */

// ---------------example--------------------
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // fprintf(stderr,"\nThis example demonstrates use of a few of the basic"
	// 	" capabilities of freesteam. It calculates the rise in temperature"
	// 	" seen in isentropic compression of a sample of steam initially at"
	// 	" 1 bar, 400 K, to a final pressure of 10 bar. It also calculates"
	// 	" the saturation temperature for steam at that final pressure.\n\n"
	// );

	// double tt = 273.15; /* in Kelvin! */
	// double pp = 1e6; /* = 1 bar */

	// fprintf(stderr,"Initial temperature = %f K, pressure = %f bar\n", tt, pp/1e5);

	// /* set a steam state of 1 bar, 400 K */
	// SteamState S = freesteam_set_pT(pp, tt);

	// Info<<"Density: "<<freesteam_rho(S)<<endl;
    // Info<<"Enthalpy: "<<freesteam_h(S)<<endl;
    // Info<<"Isobaric heat capacity: "<<freesteam_cp(S)<<endl;
    // Info<<"Isochoric heat capacity: "<<freesteam_cv(S)<<endl;
    // Info<<"Specific volume: "<<freesteam_v(S)<<endl;
    // Info<<"Internal energy: "<<freesteam_u(S)<<endl;
    // Info<<"Thermal conductivity: "<<freesteam_k(S)<<endl;
    // Info<<"Dynamic viscosity: "<<freesteam_mu(S)<<endl;
    // Info<<"Speed of sound: "<<freesteam_w(S)<<endl;
    // Info<<"Entropy: "<<freesteam_s(S)<<endl;
// ==============================================================================

#ifndef IAPWSIF97_H
#define IAPWSIF97_H

#include "steam.H"

#ifdef __cplusplus
	#define EXTERN extern "C"
#else
	#define EXTERN extern
#endif 

// pressure and temperature limitation
#define MAXT_FREESTEAM 1073.15
// if T<4C, ALPHA will be negative
#define MINT_FREESTEAM 277.15 
#define MAXP_FREESTEAM 1000E5
#define MINP_FREESTEAM 1E5
#define K2C 273.15

//CL: 
EXTERN double freesteam_p(SteamState S);
EXTERN double freesteam_T(SteamState S);
EXTERN double freesteam_rho(SteamState S);
EXTERN double freesteam_v(SteamState S);
EXTERN double freesteam_u(SteamState S);
EXTERN double freesteam_h(SteamState S);
EXTERN double freesteam_s(SteamState S);
EXTERN double freesteam_cp(SteamState S);
EXTERN double freesteam_cv(SteamState S);
EXTERN double freesteam_w(SteamState S);
EXTERN double freesteam_x(SteamState S);
EXTERN double freesteam_mu(SteamState S);
EXTERN double freesteam_k(SteamState S);

//CL: getting SteamState for two given properties e.g. pressure and temperatur
EXTERN SteamState freesteam_set_pv(double,double);
EXTERN SteamState freesteam_set_pu(double,double);
EXTERN SteamState freesteam_set_pT(double,double);
EXTERN SteamState freesteam_set_ps(double p, double s);
EXTERN SteamState freesteam_set_ph(double,double);

//CL: getting region of the SteamState
EXTERN int freesteam_region(SteamState);

//CL: transport properties
EXTERN double freesteam_mu_rhoT(double,double);
EXTERN double freesteam_k_rhoT(double,double);

//CL: Region 1 --> see region1.h (freesteam)
EXTERN double freesteam_region1_v_pT(double,double);
EXTERN double freesteam_region1_h_pT(double,double);
EXTERN double freesteam_region1_kappaT_pT(double,double);
EXTERN double freesteam_region1_alphav_pT(double,double);
EXTERN double freesteam_region1_cp_pT(double,double);
EXTERN double freesteam_region1_u_pT(double,double);
EXTERN double freesteam_region1_s_pT(double,double);
EXTERN double freesteam_region1_cv_pT(double,double);

//CL: Region 2 --> see region2.h (freesteam)
EXTERN double freesteam_region2_v_pT(double,double);
EXTERN double freesteam_region2_u_pT(double,double);
EXTERN double freesteam_region2_s_pT(double,double);
EXTERN double freesteam_region2_h_pT(double,double);
EXTERN double freesteam_region2_cp_pT(double,double);
EXTERN double freesteam_region2_cv_pT(double,double);
EXTERN double freesteam_region2_alphav_pT(double,double);
EXTERN double freesteam_region2_kappaT_pT(double,double);

//CL: Region 3 --> see region3.h (freesteam)
EXTERN double freesteam_region3_p_rhoT(double,double);
EXTERN double freesteam_region3_u_rhoT(double,double);
EXTERN double freesteam_region3_s_rhoT(double,double);
EXTERN double freesteam_region3_h_rhoT(double,double);
EXTERN double freesteam_region3_cp_rhoT(double,double);
EXTERN double freesteam_region3_cv_rhoT(double,double);
EXTERN double freesteam_region3_alphap_rhoT(double,double);
EXTERN double freesteam_region3_betap_rhoT(double,double);

//CL: Region 4 --> see region4.h (freesteam)
EXTERN double freesteam_region4_psat_T(double);
EXTERN double freesteam_region4_Tsat_p(double);
EXTERN double freesteam_region4_rhof_T(double);
EXTERN double freesteam_region4_rhog_T(double);
EXTERN double freesteam_region4_v_Tx(double,double);
EXTERN double freesteam_region4_u_Tx(double,double);
EXTERN double freesteam_region4_h_Tx(double,double);
EXTERN double freesteam_region4_s_Tx(double,double);
EXTERN double freesteam_region4_cp_Tx(double,double);
EXTERN double freesteam_region4_cv_Tx(double,double);
EXTERN double freesteam_region4_dpsatdT_T(double);

void calculateProperties_h
(
    SteamState S, 
    double &h, 
    double &rho, 
    double &mu, 
    double &kappa, 
	double &beta,
	double &cp
)
{
    int region=freesteam_region(S);
	double p,T;
    if (region==1)
    {
        p=S.R1.p;
        T=S.R1.T;
        rho=1/freesteam_region1_v_pT(S.R1.p,S.R1.T);
        h=freesteam_region1_h_pT(S.R1.p,S.R1.T);

        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        kappa=freesteam_region1_kappaT_pT(S.R1.p,S.R1.T);
        beta=freesteam_region1_alphav_pT(S.R1.p,S.R1.T);
        cp=freesteam_region1_cp_pT(S.R1.p,S.R1.T);
 
        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
    }
    else if (region==2)
    {
        p=S.R2.p;
        T=S.R2.T;
        rho=1/freesteam_region2_v_pT(S.R2.p,S.R2.T);
        h=freesteam_region2_h_pT(S.R2.p,S.R2.T);
        // x=1;

        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        kappa=freesteam_region2_kappaT_pT(S.R2.p,S.R2.T);
        beta=freesteam_region2_alphav_pT(S.R2.p,S.R2.T);
        cp=freesteam_region2_cp_pT(S.R2.p,S.R2.T);

        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
    }
    else if (region==3)
    {
        double gamma,cv;
 
        rho=S.R3.rho;
        T=S.R3.T;
        p=freesteam_region3_p_rhoT(S.R3.rho,S.R3.T);
        h=freesteam_region3_h_rhoT(S.R3.rho,S.R3.T);
        gamma=freesteam_region3_alphap_rhoT(S.R3.rho,S.R3.T);
        cp=freesteam_region3_cp_rhoT(S.R3.rho,S.R3.T);
        cv=freesteam_region3_cv_rhoT(S.R3.rho,S.R3.T);
        beta=(cp-cv)/(S.R3.T/S.R3.rho*p*gamma);
        kappa=(cp-cv)/(S.R3.T/S.R3.rho*p*p*gamma*gamma);

        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
    }
    else if (region==4)
    {
        double rhov,rhol,betav,betal,kappav,kappal,vv,vl,cpl,cpv,hl,hv,cp;
        double dvldp,dvvdp,dhldp,dhvdp;
        double dpdT,dvdh,dvdp,dxdp;

        SteamState Sl,Sv;
      
        // x=S.R4.x;
        T=S.R4.T;
        rho=1/freesteam_region4_v_Tx(S.R4.T,S.R4.x);
        h=freesteam_region4_h_Tx(S.R4.T,S.R4.x);
        p=freesteam_region4_psat_T(S.R4.T);
        cp=freesteam_region4_cp_Tx(S.R4.T,S.R4.x);
        //CL: Getting density on the vapour and liquid lines
        rhov=freesteam_region4_rhog_T(S.R4.T);
        rhol=freesteam_region4_rhof_T(S.R4.T);
        vv=1/rhov;
        vl=1/rhol;
        //CL: getting derivatives --> this is a bit tricky inside the vapor dome
        dpdT=freesteam_region4_dpsatdT_T(S.R4.T);
	    // getting the states outside the vapour dome 
        Sl=freesteam_set_pv(p,vl-0.0000001);  //inside region 1
        Sv=freesteam_set_pv(p,vv+0.0000001);  //inside region 2
  
        kappal=freesteam_region1_kappaT_pT(Sl.R1.p,Sl.R1.T);
        kappav=freesteam_region2_kappaT_pT(Sv.R2.p,Sv.R2.T);

        betal=freesteam_region1_alphav_pT(Sl.R1.p,Sl.R1.T);
        betav=freesteam_region2_alphav_pT(Sv.R2.p,Sv.R2.T);

        cpl=freesteam_region1_cp_pT(Sl.R1.p,Sl.R1.T);
        cpv=freesteam_region2_cp_pT(Sv.R2.p,Sv.R2.T);

        hl=freesteam_region1_h_pT(Sl.R1.p,Sl.R1.T);
        hv=freesteam_region2_h_pT(Sv.R2.p,Sv.R2.T);
        //calculation derviatives on liquid and vapour line
        dvldp=betal*vl/dpdT-kappal*vl;
        dvvdp=betav*vv/dpdT-kappav*vv;
        dhldp=vl*(1-betal*Sl.R1.T)+cpl/dpdT;
        dhvdp=vv*(1-betav*Sv.R2.T)+cpl/dpdT;
        dxdp=-dhldp/(hv-hl)
                 +(h-hl)/((hv-hl)*(hv-hl))
                     *(dhvdp-dhldp);
        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
    }
    else
    {
        Info<<"IAPWS-IF97 error, outside the regions 1-4"<<endl;
    }
}

void calculateRho_h(SteamState S, double &rho)
{
    int region=freesteam_region(S);
    if (region==1)
    {
        rho=1/freesteam_region1_v_pT(S.R1.p,S.R1.T);
    }
    else if (region==2)
    {
        rho=1/freesteam_region2_v_pT(S.R2.p,S.R2.T);
    }
    else if (region==3)
    {
        rho=S.R3.rho;
    }
    else if (region==4)
    {
        rho=1/freesteam_region4_v_Tx(S.R4.T,S.R4.x);
    }
    else
    {
        Info<<"IAPWS-IF97 error, outside the regions 1-4"<<endl;
    }
}

void calculateDlnRho_DlnT(double p,double T, double& dLnRho_dLnT,double dT=1)
{
    double T1=T-dT/2.0, T2=T+dT/2.0;
    SteamState S1 = freesteam_set_pT(p, T1);
    SteamState S2 = freesteam_set_pT(p, T2);
    double rho1, rho2;
    calculateRho_h(S1, rho1);
    calculateRho_h(S2, rho2);
    double lnRho1=std::log(rho1);
    double lnRho2=std::log(rho2);
    double lnT1=std::log(T1);
    double lnT2=std::log(T2);
    dLnRho_dLnT=(lnRho2-lnRho1)/(lnT2-lnT1);
}
#endif