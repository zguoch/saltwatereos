#!python
# -*-coding:utf-8-*-
# Calculate phase diagram and thermal dynamic properties of NaCl-H2O system in 2D space
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/01/05, GEOMAR
# ===============================================================

import sys
import argparse
import os
#===============================================================
import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
# config font
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 12
mpl.rcParams['mathtext.fontset'] = 'cm'
sec_year=86400*365
cm2inch=1 #0.393701
dpi=600

from pyswEOS import H2ONaCl
from pyswEOS import H2O
water=H2O.cH2O()
# sw=H2ONaCl.cH2ONaCl()
fmt_figs=['pdf']
dpi=400
def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Calculate phase diagram and thermal dynamic properties of NaCl-H2O system in 2D space'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/01/05, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))

def XT_P0(P0=300, Xmin=1E-5, Xmax=0.999999, Tmin=1, Tmax=1000, nX=100, nT=100):
    """Calculate phase diagram and properties in X-T space 
    with given P

    Args:
        P0 (float, optional): Pressure [bar]. Defaults to 100.
        Xmin (float, optional): Min salinity [wt% NaCl]. Defaults to 0.
        Xmax (float, optional): Max salinity [wt% NaCl]. Defaults to 1.
        Tmin (float, optional): Min temperature [deg. C]. Defaults to 0.
        Tmax (float, optional): Max temperature [deg. C]. Defaults to 500.
    """
    Hmax = 4 # 4 MJ/kg
    sw=H2ONaCl.cH2ONaCl()
    T=np.linspace(Tmin, Tmax, nT)
    X=np.linspace(Xmin, Xmax, nX)
    P_vlh=np.zeros_like(T)
    T_vlh=np.array(sw.T_VaporLiquidHaliteCoexist(P0))
    T_boil_water=water.T_Boiling(P0)
    X_LH=np.zeros_like(T)
    # 1. NaCl liquidus, only valid when P > P_{V+L+H}
    for i in range(0,len(T)):
        P_vlh[i] = sw.P_VaporLiquidHaliteCoexist(T[i]);
    if(P0>389):
        T_vlh=[]
    # plot
    # fig=plt.figure()
    fig,axes=plt.subplots(1,2,gridspec_kw={"width_ratios":[1,1],'wspace':0.1},figsize=(12,5))
    ax_T,ax_H=axes[0],axes[1]
    ax_H.yaxis.set_ticks_position('right')
    ax_H.yaxis.set_label_position('right')
    # 1. V+L+H
    # region contourf in PTX space
    # XX,TT=np.meshgrid(X,T)
    # Regions_T=np.zeros_like(XX)
    # for i in range(0,len(T)):
    #     for j in range(0,len(X)):
    #         prop=sw.prop_pTX(P0*1E5, T[i]+273.15, X[j])
    #         Regions_T[i][j]=prop.Region
    # ax_T.contourf(XX,TT, Regions_T, levels=50)
    # region contourf in PHX space
    # H=np.linspace(1E-3,4,100)
    # XX,HH=np.meshgrid(X,H)
    # Regions=np.zeros_like(XX)
    # for i in range(0,len(H)):
    #     for j in range(0,len(X)):
    #         prop=sw.prop_pHX(P0*1E5, H[i]*1E6, X[j])
    #         Regions[i][j]=prop.Region
    # ax_H.contourf(XX,HH, Regions, levels=50)
    # 2. L+H
    if(len(T_vlh)==2):
        T_vlh[0]=T_vlh[0]-1 # remove computation error, ensure it doesn't go cross LVH region
        T_vlh[1]=T_vlh[1]+1
        # segment 1
        # PTX space
        T_LH1=np.linspace(Tmin, T_vlh[0], 50)
        X_LH1=sw.Mol2Wt(sw.X_HaliteLiquidus(T_LH1, T_LH1*0+P0))
        l_LH,=ax_T.plot(X_LH1,T_LH1, color='g')
        # LVH point
        p_LVH=ax_T.scatter(X_LH1[-1],T_LH1[-1],marker='o', zorder=10, fc='orange',ec='blue',label='LVH')

        H_LH1=np.zeros_like(T_LH1)
        for i in range(0,len(T_LH1)):
            prop=sw.prop_pTX(P0*1E5, T_LH1[i]+273.15, X_LH1[i])
            H_LH1[i]=prop.H
        # PHX space
        ax_H.plot(X_LH1,H_LH1/1E6,color=l_LH.get_color())
        prop=sw.prop_pHX(P0*1E5, H_LH1[-1], X_LH1[-1])
        # LVH region
        dH=1E4 #0.01 MJ/kg
        while(not prop.Region==5):
            prop=sw.prop_pHX(P0*1E5, H_LH1[-1]+dH, X_LH1[-1])
            dH+=1E4
            if(np.abs(dH)>1E5):
                break
        ax_H.fill([prop.X_l, prop.X_v, 1, prop.X_l], 
                np.array([prop.H_l, prop.H_v, prop.H_h, prop.H_l])/1E6,
                fc=p_LVH._facecolors[0], ec=p_LVH._edgecolors[0], label='LVH')
        H_VH1=np.array([prop.H_v, prop.H_h])
        # segment 2
        T_LH2=np.linspace(T_vlh[1],Tmax, 50)
        X_LH2=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(T_LH2, T_LH2*0+P0)))
        ax_T.plot(X_LH2,T_LH2, color='g')
        # LVH point
        ax_T.plot(X_LH2[0],T_LH2[0],marker='o', zorder=10, mfc=p_LVH._facecolors[0],mec=p_LVH._edgecolors[0])
        H_LH2=np.zeros_like(T_LH2)
        for i in range(0,len(T_LH2)):
            prop=sw.prop_pTX(P0*1E5, T_LH2[i]+273.15, X_LH2[i])
            H_LH2[i]=prop.H
        # PHX space
        ax_H.plot(X_LH2,H_LH2/1E6,color=l_LH.get_color())
        prop=sw.prop_pHX(P0*1E5, H_LH2[0], X_LH2[0])
        # LVH region
        dX=0.005
        while(not prop.Region==5):
            prop=sw.prop_pHX(P0*1E5, H_LH2[0], X_LH2[0]+dX)
            # dH-=dH
            dX-=0.001
            if(np.abs(dX)>0.005):
                break
        ax_H.fill([prop.X_l, prop.X_v, 1, prop.X_l], 
                np.array([prop.H_l, prop.H_v, prop.H_h, prop.H_l])/1E6,
                fc=p_LVH._facecolors[0], ec=p_LVH._edgecolors[0])
        H_VH2=np.array([prop.H_v, prop.H_h])
        
        # V + H
        fill_VH=ax_T.axhspan(T_vlh[0],T_vlh[1],fc='lightgray',ec=p_LVH._edgecolors[0],alpha=0.5,label='V+H')
        X_VH=np.array([0, 1])
        ax_H.fill_between(X_VH, H_VH1/1E6, H_VH2/1E6, ec=p_LVH._edgecolors[0], fc=fill_VH._facecolor,label='V+H')
    else:
        X_LH=sw.Mol2Wt(sw.X_HaliteLiquidus(T, T*0+P0))
        l_LH,=ax_T.plot(X_LH,T,color='g')
        H_LH=np.zeros_like(T)
        for i in range(0,len(T)):
            prop=sw.prop_pTX(P0*1E5, T[i]+273.15, X_LH[i])
            H_LH[i]=prop.H
        ax_H.plot(X_LH,H_LH/1E6,color=l_LH.get_color())
    # 3. L+V
    T_crit_P0, X_crit_P0=sw.T_X_Critical(P0)
    prop_crit=sw.prop_pTX(P0*1E5, T_crit_P0+273.15, sw.Mol2Wt(X_crit_P0))
    if(len(T_vlh)==2):
        # segment 1
        T_LV_L1=np.linspace(T_crit_P0,T_vlh[0],100)
        X_LV_L1=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(T_LV_L1,T_LV_L1*0+P0)))
        T_LV_V=np.linspace(T_crit_P0,Tmax,100)
        X_LV_V=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_VaporBranch(T_LV_V,T_LV_V*0+P0)))
        l_LV_L,=ax_T.plot(X_LV_L1, T_LV_L1, color='k', label='L+V L')
        l_LV_V,=ax_T.plot(X_LV_V, T_LV_V, color='r', label='L+V V')
        p_Critical,=ax_T.plot(sw.Mol2Wt(X_crit_P0), T_crit_P0, marker='o', zorder=10,  mfc='r', mec='orange')
        p_Boil_water,=ax_T.plot(0,T_boil_water,marker='o',zorder=10, mfc='r',mec='blue')
        H_LV_L1=np.zeros_like(T_LV_L1)
        H_LV_V=np.zeros_like(T_LV_V)
        for i in range(0,len(T_LV_L1)):
            prop=sw.prop_pTX(P0*1E5, T_LV_L1[i]+273.15, X_LV_L1[i])
            H_LV_L1[i]=prop.H
        for i in range(0,len(T_LV_V)):
            prop=sw.prop_pTX(P0*1E5, T_LV_V[i]+273.15, X_LV_V[i])
            H_LV_V[i]=prop.H
        ax_H.plot(X_LV_L1,H_LV_L1/1E6,color=l_LV_L.get_color())
        ax_H.plot(X_LV_V,H_LV_V/1E6,color=l_LV_V.get_color())
        ax_H.plot(sw.Mol2Wt(X_crit_P0), prop_crit.H/1E6, marker='o', zorder=10, mfc=p_Critical._markerfacecolor,mec=p_Critical._markeredgecolor)
        # prop=sw.prop_pTX(P0*1E5, T_boil_water+273.15, 0)
        # ax_H.plot(0, prop.H_l/1E6, marker='o', zorder=10, mfc=p_Boil_water._markerfacecolor,mec=p_Boil_water._markeredgecolor)
        # segment 2
        T_LV_L2=np.linspace(T_vlh[1],Tmax,100)
        X_LV_L2=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(T_LV_L2,T_LV_L2*0+P0)))
        ax_T.plot(X_LV_L2, T_LV_L2, color=l_LV_L.get_color())
        H_LV_L2=np.zeros_like(T_LV_L2)
        for i in range(0,len(T_LV_L2)):
            prop=sw.prop_pTX(P0*1E5, T_LV_L2[i]+273.15, X_LV_L2[i])
            H_LV_L2[i]=prop.H
        ax_H.plot(X_LV_L2,H_LV_L2/1E6,color=l_LV_L.get_color())
        
    else:
        T_LV1=np.linspace(T_crit_P0,Tmax,100)
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(T_LV1,T_LV1*0+P0)))
        X_LV_V=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_VaporBranch(T_LV1,T_LV1*0+P0)))
        l_LV_L,=ax_T.plot(X_LV_L, T_LV1, color='k')
        l_LV_V,=ax_T.plot(X_LV_V, T_LV1, color='red')
        p_Critical,=ax_T.plot(sw.Mol2Wt(X_crit_P0), T_crit_P0, marker='o', zorder=10,  mfc='r', mec='orange')
        
        H_LV_L=np.zeros_like(T_LV1)
        H_LV_V=np.zeros_like(T_LV1)
        for i in range(0,len(T_LV1)):
            prop=sw.prop_pTX(P0*1E5, T_LV1[i]+273.15, X_LV_L[i])
            H_LV_L[i]=prop.H
            prop=sw.prop_pTX(P0*1E5, T_LV1[i]+273.15, X_LV_V[i])
            H_LV_V[i]=prop.H
        ax_H.plot(X_LV_L,H_LV_L/1E6,color=l_LV_L.get_color())
        ax_H.plot(X_LV_V,H_LV_V/1E6,color=l_LV_V.get_color())
        ax_H.plot(sw.Mol2Wt(X_crit_P0), prop_crit.H/1E6, marker='o', zorder=10, mfc=p_Critical._markerfacecolor,mec=p_Critical._markeredgecolor)
    # upper and lower boundary of enthalpy
    H_boundary_Tmax,H_boundary_Tmin=np.zeros_like(X),np.zeros_like(X)
    for i in range(0,len(X)):
        prop=sw.prop_pTX(P0*1E5, Tmax+273.15, X[i])
        H_boundary_Tmax[i]=prop.H
        prop=sw.prop_pTX(P0*1E5, Tmin+273.15, X[i])
        H_boundary_Tmin[i]=prop.H
    ax_H.plot(X,H_boundary_Tmax/1E6, 'k')
    fill_invalid_H=ax_H.fill_between(X, H_boundary_Tmax/1E6, Hmax, color='k', alpha=0.5, label='Invalid region')
    ax_H.fill_between(X, H_boundary_Tmin/1E6, 0, color='k', alpha=0.5)
    
    ax_H.legend()
    ax_T.legend()
    # axis
    for ax in axes:
        ax.set_xlim(Xmin, Xmax)
        ax.set_xlabel('Salinity (wt% NaCl)')
        ticks=ax.get_xticks()
        ax.xaxis.set_ticks(ticks)
        labels=[]
        for tick in ticks:
            labels.append(str('%.0f'%(tick*100)))
        ax.xaxis.set_ticklabels(labels)
    ax_T.set_ylim(Tmin,Tmax)
    ax_T.set_ylabel('Temperature ($^{\circ}$C)')
    ax_H.set_ylabel('Specific enthalpy (MJ/kg)')
    ax_H.set_ylim(0,Hmax)
    # ax_H.grid(axis='both',which='both')

    plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.98)
    for fmt in fmt_figs:
        plt.savefig('PhaseDiagram_XT_P0.%s'%(fmt), dpi=dpi)
    # plt.show()

def main(argv):
    # argc=len(argv)
    # usage(argv)
    # exit(0)

    XT_P0()

if __name__ == '__main__':
    sys.exit(main(sys.argv))