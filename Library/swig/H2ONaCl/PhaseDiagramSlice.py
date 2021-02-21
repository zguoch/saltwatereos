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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
# config font
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 12
mpl.rcParams['mathtext.fontset'] = 'cm'
sec_year=86400*365
cm2inch=1 #0.393701
dpi=600

from pyswEOS import H2ONaCl
from pyswEOS import H2O
from pyswEOS import NaCl
water=H2O.cH2O()
sw=H2ONaCl.cH2ONaCl()
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
def H_PTX(P_bar, T_C, X_wt):
    H=np.zeros_like(P_bar)
    for j in range(0,len(H)):
        prop=sw.prop_pTX(P_bar[j]*1E5, T_C[j]+273.15, X_wt[j])
        H[j]=prop.H
    return H
def XT_P0(P0=300, Xmin=1E-5, Xmax=0.999999, Tmin=1, Tmax=1000, nX=100, nT=100, 
    fc_VLH=(1,1,0,0.8), ec_VLH='b',fc_VH=(0.5,0,0.5,0.5), lc_LH='b',lc_LV_L='b',
    lc_LV_V='r',fc_LH='lightblue',fc_LV='lightgreen'):
    """Calculate phase diagram and properties in X-T space 
    with given P

    Args:
        P0 (float, optional): Pressure [bar]. Defaults to 100.
        Xmin (float, optional): Min salinity [wt% NaCl]. Defaults to 0.
        Xmax (float, optional): Max salinity [wt% NaCl]. Defaults to 1.
        Tmin (float, optional): Min temperature [deg. C]. Defaults to 0.
        Tmax (float, optional): Max temperature [deg. C]. Defaults to 500.
    """
    Hmax = 4E6 # 4 MJ/kg
    zorder_p=10
    sw=H2ONaCl.cH2ONaCl()
    T=np.linspace(Tmin, Tmax, nT)
    X=np.linspace(Xmin, Xmax, nX)
    # 1.0 upper and lower boundary of enthalpy
    H_boundary_Tmax,H_boundary_Tmin=np.zeros_like(X),np.zeros_like(X)
    for i in range(0,len(X)):
        prop=sw.prop_pTX(P0*1E5, Tmax+273.15, X[i])
        H_boundary_Tmax[i]=prop.H
        prop=sw.prop_pTX(P0*1E5, Tmin+273.15, X[i])
        H_boundary_Tmin[i]=prop.H
    # 1.1 water boiling temperature 
    T_boil_water=water.T_Boiling(P0)
    # 1.2 V+L+H enthalpy and salinity when 1<P<390 bar
    HminHmaxXminXmax=sw.HX_VaporLiquidHaliteCoexist(P0) 
    HX=HminHmaxXminXmax
    # 1.3 critical point
    T_crit_P0, X_crit_P0=sw.T_X_Critical(P0)
    prop_crit=sw.prop_pTX(P0*1E5, T_crit_P0+273.15, sw.Mol2Wt(X_crit_P0))
    # 1.4 Phase boundary of L+H
    T_LH, T_LV_V = [np.linspace(Tmin, Tmax, 100)], np.linspace(T_crit_P0, Tmax, 100)
    T_LV_L=[T_LV_V]
    # plot
    # fig=plt.figure()
    fig,axes=plt.subplots(1,2,gridspec_kw={"width_ratios":[1,1],'wspace':0.1},figsize=(12,5))
    ax_T,ax_H=axes[0],axes[1]
    ax_H.yaxis.set_ticks_position('right')
    ax_H.yaxis.set_label_position('right')
    # 2.
    if(len(HX)==4):
        prop1=sw.prop_pHX(P0*1E5, HX[0], HX[2])
        prop2=sw.prop_pHX(P0*1E5, HX[1], HX[3])
        # 2.1 PTX space
        # L+V+H
        p_LVH=ax_T.scatter(HX[2:], np.array([prop1.T,prop2.T]), fc=fc_VLH, ec=ec_VLH,zorder=zorder_p,label='V+L+H')
        # VH
        fill_VH=ax_T.axhspan(prop1.T, prop2.T, fc=fc_VH, ec=ec_VLH, label='V+H')
        # 2.2 PHX space
        # VLH
        ax_H.fill_between([Xmin,Xmax], [prop1.H_v, prop1.H_h], [prop2.H_v, prop2.H_h], fc=fc_VH, ec=ec_VLH, label='V+H')
        ax_H.fill([Xmin, HX[2], Xmax, Xmin],[prop1.H_v, prop1.H_l, prop1.H_h, prop1.H_v], fc=fc_VLH,ec=ec_VLH, label='L+V+H')
        ax_H.fill([Xmin, HX[3], Xmax, Xmin],[prop2.H_v, prop2.H_l, prop2.H_h, prop2.H_v], fc=fc_VLH,ec=ec_VLH)
        
        # LH
        T_LH=[np.linspace(Tmin, prop1.T-1, 100), np.linspace(prop2.T+1, Tmax, 100)]
        # first part
        TT=T_LH[0]
        X_LH=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(TT, TT*0+P0)))
        ax_T.plot(X_LH, TT, color=lc_LH)
        ax_T.fill_betweenx(TT, X_LH, Tmax, fc=fc_LH, ec='None', alpha=0.8,label='L+H')
        H_LH=H_PTX(P0+TT*0, TT, X_LH)
        ax_H.plot(X_LH, H_LH, color=lc_LH)
        X_LH_bottom=np.linspace(Xmax, X_LH.min(), 10)
        H_LH_bottom=H_PTX(P0+X_LH_bottom*0, Tmin+X_LH_bottom*0, X_LH_bottom)
        ax_H.fill(np.append(np.append(X_LH_bottom, X_LH),(Xmax,Xmax)), np.append(np.append(H_LH_bottom, H_LH),(prop1.H_h, H_LH_bottom[0])),
                fc=fc_LH, ec='None', alpha=0.8,label='L+H')
        # second part
        TT=T_LH[1]
        X_LH=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(TT, TT*0+P0)))
        ax_T.plot(X_LH, TT, color=lc_LH)
        ax_T.fill_betweenx(TT, X_LH, Tmax, fc=fc_LH, ec='None', alpha=0.8)
        H_LH=H_PTX(P0+TT*0, TT, X_LH)
        ax_H.plot(X_LH, H_LH, color=lc_LH)
        ax_H.fill(np.append(X_LH,(Xmax,X_LH[0])), np.append(H_LH,(prop2.H_h, H_LH[0])),
                fc=fc_LH, ec='None', alpha=0.8)
        # LV_L
        # first part
        T_LV_L=[np.linspace(T_crit_P0, prop1.T-1, 100), np.linspace(prop2.T+1, Tmax, 100)]
        TT=T_LV_L[0]
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(TT,TT*0+P0)))
        ax_T.plot(X_LV_L, TT, color=lc_LV_L)
        ax_T.fill_between(X_LV_L, TT, prop1.T, fc=fc_LV, ec='None', alpha=0.8, label='L+V')
        X_LV_V=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_VaporBranch(T_LV_V,T_LV_V*0+P0)))
        ax_T.plot(X_LV_V, T_LV_V, color=lc_LV_V)
        ax_T.fill_between(X_LV_V, T_LV_V, prop1.T, fc=fc_LV, ec='None', alpha=0.8)
        H_LV_L=H_PTX(P0+TT*0, TT, X_LV_L)
        ax_H.plot(X_LV_L, H_LV_L, color=lc_LV_L)
        H_LV_V=H_PTX(P0+TT*0, T_LV_V, X_LV_V)
        ax_H.plot(X_LV_V, H_LV_V, color=lc_LV_V)
        ax_H.fill(np.append(np.append(np.flip(X_LV_V),X_LV_L),Xmin), 
                np.append(np.append(np.flip(H_LV_V),H_LV_L),prop1.H_v), 
                fc=fc_LV, ec='None', alpha=0.8, label='L+V')
        # print(H_LV_L,H_LV_V)
        # second part
        TT=T_LV_L[1]
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(TT,TT*0+P0)))
        ax_T.plot(X_LV_L, TT, color=lc_LV_L)
        ax_T.fill_betweenx(TT, X_LV_L, Xmin, fc=fc_LV, ec='None', alpha=0.8)
        H_LV_L=H_PTX(P0+TT*0, TT, X_LV_L)
        ax_H.plot(X_LV_L, H_LV_L, color=lc_LV_L)
        ax_H.fill_between(np.insert(X_LV_L,0,Xmin), np.insert(H_LV_L,0,prop2.H_v), Hmax, fc=fc_LV, ec='None', alpha=0.8)

    else:
        TT=T_LH[0]
        X_LH=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(TT, TT*0+P0)))
        ax_T.plot(X_LH, TT, color=lc_LH)
        ax_T.fill_betweenx(TT, X_LH, Tmax, fc=fc_LH, ec='None', alpha=0.8,label='L+H')
        H_LH=H_PTX(P0+TT*0, TT, X_LH)
        ax_H.plot(X_LH, H_LH, color=lc_LH)
        H_LH_bottom=H_PTX(P0+X_LH*0, Tmin+X_LH*0, X_LH)
        ax_H.fill_between(X_LH, H_LH, H_LH_bottom,fc=fc_LH, ec='None', alpha=0.8,label='L+H')
        # LV
        TT=T_LV_L[0]
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(TT,TT*0+P0)))
        X_LV_L[0]=sw.Mol2Wt(X_crit_P0)
        ax_T.plot(X_LV_L, TT, color=lc_LV_L)
        ax_T.fill_between(X_LV_L, TT, Tmax, fc=fc_LV, ec='None', alpha=0.8, label='L+V')
        X_LV_V=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_VaporBranch(T_LV_V,T_LV_V*0+P0)))
        X_LV_V[0]=sw.Mol2Wt(X_crit_P0)
        ax_T.plot(X_LV_V, T_LV_V, color=lc_LV_V)
        ax_T.fill_between(X_LV_V, T_LV_V, Tmax, fc=fc_LV, ec='None', alpha=0.8)
        H_LV_L=H_PTX(P0+TT*0, TT, X_LV_L)
        ax_H.plot(X_LV_L, H_LV_L, color=lc_LV_L)
        ax_H.fill_between(X_LV_L, H_LV_L, Hmax, fc=fc_LV, ec='None', alpha=0.8, label='L+V')
        H_LV_V=H_PTX(P0+TT*0, T_LV_V, X_LV_V)
        ax_H.plot(X_LV_V, H_LV_V, color=lc_LV_V)
        ax_H.fill_between(X_LV_V, H_LV_V, Hmax, fc=fc_LV, ec='None', alpha=0.8)
    # critical point
    ax_T.scatter(sw.Mol2Wt(X_crit_P0), T_crit_P0, fc='r',ec='orange', zorder=zorder_p, label='Critical point')
    ax_H.scatter(sw.Mol2Wt(X_crit_P0), prop_crit.H, fc='r',ec='orange', zorder=zorder_p, label='Critical point')
    # valid region
    fill_invalid=ax_H.fill_between(X, H_boundary_Tmax, Hmax, fc='lightgray', ec='k', label='Invalid region')
    ax_H.fill_between(X, H_boundary_Tmin, 0, fc='lightgray', ec='k')
    ax_H.axhspan(1,2,fc='w',ec='b',label='L')
    ax_T.axhspan(1,2,fc='w',ec='b',label='L')
    ax_T.legend(ncol=2)
    ax_H.legend(ncol=2)
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
    # enthalpy axis
    ticks=ax_H.get_yticks()
    ax_H.yaxis.set_ticks(ticks)
    labels=[]
    for tick in ticks:
        labels.append(str('%.1f'%(tick/1E6)))
    ax_H.yaxis.set_ticklabels(labels)
    ax_H.set_ylim(0,Hmax)
    # ax_H.grid(axis='both',which='both')

    plt.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.98)
    for fmt in fmt_figs:
        plt.savefig('PhaseDiagram_XT_P0_%.0fbar.%s'%(P0,fmt), dpi=dpi)
    # plt.show()
line_intersection = lambda p1, p2, x0: (x0-p1[0])/(p2[0]-p1[0])*(p2[1]-p1[1])+p1[1]
def XT_P0_logLinearScale(P0=150, X0=0.2, H0=1.0E6, Xmin=0.5E-4, Xmax=0.999999, Tmin=5, Tmax=1000, nX=100, nT=100, 
    fc_VLH=(1,1,0,0.8), ec_VLH='k',fc_VH=(0.5,0,0.5,0.5), lc_LH='k',lc_LV_L='k',
    lc_LV_V='r',fc_LH='lightblue',fc_LV='lightgreen', savefig=True,lw=0.5,alpha=1,X_max_log=0.1, Hmax=4.2E6):
    """Calculate phase diagram and properties in X-T space 
    with given P

    Args:
        P0 (float, optional): Pressure [bar]. Defaults to 100.
        X0 (float, optional): bulk composition [wt% NaCl]
        H0 (float, optional): phase state with constant bulk enthalpy [J/kg]
        Xmin (float, optional): Min salinity [wt% NaCl]. Defaults to 0.
        Xmax (float, optional): Max salinity [wt% NaCl]. Defaults to 1.
        Tmin (float, optional): Min temperature [deg. C]. Defaults to 0.
        Tmax (float, optional): Max temperature [deg. C]. Defaults to 500.
    """
    # Hmax = 4E6 # 4 MJ/kg
    zorder_p=10
    sw=H2ONaCl.cH2ONaCl()
    X_max_log=X_max_log
    T=np.linspace(Tmin, Tmax, nT)
    X_log=np.linspace(np.log10(Xmin), np.log10(X_max_log*1.1), 20)
    X_log=10**X_log
    X_linear=np.linspace(X_max_log, Xmax, nX-len(X_log))
    # X=np.linspace(Xmin, Xmax, nX)
    X=np.append(X_log, X_linear)
    # 1.0 upper and lower boundary of enthalpy
    H_boundary_Tmax,H_boundary_Tmin=np.zeros_like(X),np.zeros_like(X)
    for i in range(0,len(X)):
        prop=sw.prop_pTX(P0*1E5, Tmax+273.15, X[i])
        H_boundary_Tmax[i]=prop.H
        prop=sw.prop_pTX(P0*1E5, Tmin+273.15, X[i])
        H_boundary_Tmin[i]=prop.H
    # 1.1 water boiling temperature 
    T_boil_water=water.T_Boiling(P0)
    # 1.2 V+L+H enthalpy and salinity when 1<P<390 bar
    HminHmaxXminXmax=sw.HX_VaporLiquidHaliteCoexist(P0) 
    HX=HminHmaxXminXmax
    # 1.3 critical point
    T_crit_P0, X_crit_P0=sw.T_X_Critical(P0)
    prop_crit=sw.prop_pTX(P0*1E5, T_crit_P0+273.15, sw.Mol2Wt(X_crit_P0))
    # 1.4 Phase boundary of L+H
    T_LH, T_LV_V = [np.linspace(Tmin, Tmax, 100)], np.linspace(T_crit_P0, Tmax, 100)
    T_LV_L=[T_LV_V]
    # plot
    # fig=plt.figure()
    fig,axes=plt.subplots(1,2,gridspec_kw={"width_ratios":[1,1],'wspace':0.05},figsize=(14,6))
    ax_T,ax_H=axes[0],axes[1]
    ax_H.yaxis.set_ticks_position('right')
    ax_H.yaxis.set_label_position('right')
    # hybrid scale 
    divider = make_axes_locatable(ax_H)
    ax_H_log = divider.append_axes("left", size=2.0, pad=0)
    ax_H_log.set_xscale('log')
    ax_H_log.set_xlim((Xmin, X_max_log))
    ax_H_log.spines['right'].set_visible(False)
    ax_H.spines['left'].set_visible(False)
    ax_H_log.set_yticks([])
    ax_H_log.axvspan(0,1,fc='lightgray',alpha=0.4,zorder=10)
    # ax_H_log.text(0.5,0.01, 'Log scale',ha='center',va='bottom',color='b', transform=ax_H_log.transAxes)
    ax_H_log.text(0.5,0.1, 'H=%.2f MJ/kg'%(H0/1E6),ha='center',va='bottom',color='k',fontweight='bold', transform=ax_H_log.transAxes,zorder=10)
    ax_H_log.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=20))
    ax_H_log.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8),numticks=20))
    ax_H_log.tick_params(axis='x',which='both', colors='blue')
    ax_H_log.spines['bottom'].set_color('b')
    # 2.
    if(len(HX)==4):
        prop1=sw.prop_pHX(P0*1E5, HX[0], HX[2])
        prop2=sw.prop_pHX(P0*1E5, HX[1], HX[3])
        # 2.1 PTX space
        # L+V+H
        p_LVH=ax_T.scatter(HX[2:], np.array([prop1.T,prop2.T]), fc=fc_VLH, ec=ec_VLH,zorder=zorder_p,label='V+L+H')
        # VH
        fill_VH=ax_T.axhspan(prop1.T, prop2.T, fc=fc_VH, ec=ec_VLH, label='V+H')
        # 2.2 PHX space
        # VLH
        ax_H.fill_between([Xmin,Xmax], [prop1.H_v, prop1.H_h], [prop2.H_v, prop2.H_h], fc=fc_VH, ec=ec_VLH, lw=lw,label='V+H')
        ax_H.fill([Xmin, HX[2], Xmax, Xmin],[prop1.H_v, prop1.H_l, prop1.H_h, prop1.H_v], fc=fc_VLH,ec=ec_VLH, lw=lw, label='L+V+H')
        ax_H.fill([Xmin, HX[3], Xmax, Xmin],[prop2.H_v, prop2.H_l, prop2.H_h, prop2.H_v], fc=fc_VLH,ec=ec_VLH, lw=lw)
        # log scale branch
        # VH
        H_inter1_VH=line_intersection((Xmin, prop1.H_v),(Xmax,prop1.H_h),X_max_log)
        H_inter2_VH=line_intersection((Xmin, prop2.H_v),(Xmax,prop2.H_h),X_max_log)
        H_log_VH1=line_intersection((Xmin, prop1.H_v),(Xmax,prop1.H_h),X_log)
        H_log_VH2=line_intersection((Xmin, prop2.H_v),(Xmax,prop2.H_h),X_log)
        # ax_H_log.plot(X_log,H_log_VH1,'.')
        ax_H_log.fill_between(X_log, H_log_VH1, H_log_VH2, fc=fc_VH, ec=ec_VLH, lw=lw, label='V+H')
        # VLH
        H_log_VLH1=line_intersection((Xmin, prop1.H_v),(HX[2],prop1.H_l),X_log)
        ax_H_log.fill_between(X_log, H_log_VH1, H_log_VLH1, fc=fc_VLH,ec=ec_VLH, lw=lw, label='L+V+H')
        H_log_VLH2=line_intersection((Xmin, prop2.H_v),(HX[3],prop2.H_l),X_log)
        ax_H_log.fill_between(X_log, H_log_VH2, H_log_VLH2, fc=fc_VLH,ec=ec_VLH, lw=lw)
        
        # LH
        T_LH=[np.linspace(Tmin, prop1.T-1, 100), np.linspace(prop2.T+1, Tmax, 100)]
        # first part
        TT=T_LH[0]
        X_LH=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(TT, TT*0+P0)))
        ax_T.plot(X_LH, TT, color=lc_LH)
        ax_T.fill_betweenx(TT, X_LH, Tmax, fc=fc_LH, ec='None', alpha=alpha,label='L+H')
        H_LH=H_PTX(P0+TT*0, TT, X_LH)
        ax_H.plot(X_LH, H_LH, color=lc_LH, lw=lw)
        X_LH_bottom=np.linspace(Xmax, X_LH.min(), 10)
        H_LH_bottom=H_PTX(P0+X_LH_bottom*0, Tmin+X_LH_bottom*0, X_LH_bottom)
        ax_H.fill(np.append(np.append(X_LH_bottom, X_LH),(Xmax,Xmax)), np.append(np.append(H_LH_bottom, H_LH),(prop1.H_h, H_LH_bottom[0])),
                fc=fc_LH, ec='None', lw=lw, alpha=alpha,label='L+H')
        # second part
        TT=T_LH[1]
        X_LH=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(TT, TT*0+P0)))
        ax_T.plot(X_LH, TT, color=lc_LH)
        ax_T.fill_betweenx(TT, X_LH, Tmax, fc=fc_LH, ec='None', alpha=alpha)
        H_LH=H_PTX(P0+TT*0, TT, X_LH)
        ax_H.plot(X_LH, H_LH, color=lc_LH, lw=lw)
        ax_H.fill(np.append(X_LH,(Xmax,X_LH[0])), np.append(H_LH,(prop2.H_h, H_LH[0])),
                fc=fc_LH, ec='None', lw=lw, alpha=alpha)
        # LV_L
        # first part
        T_LV_L=[np.linspace(T_crit_P0, prop1.T-1, 100), np.linspace(prop2.T+1, Tmax, 100)]
        TT=T_LV_L[0]
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(TT,TT*0+P0)))
        ax_T.plot(X_LV_L, TT, color=lc_LV_L)
        ax_T.fill_between(X_LV_L, TT, prop1.T, fc=fc_LV, ec='None', alpha=alpha, label='L+V')
        X_LV_V=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_VaporBranch(T_LV_V,T_LV_V*0+P0)))
        ax_T.plot(X_LV_V, T_LV_V, color=lc_LV_V)
        ax_T.fill_between(X_LV_V, T_LV_V, prop1.T, fc=fc_LV, ec='None', alpha=alpha)
        H_LV_L=H_PTX(P0+TT*0, TT, X_LV_L)
        ax_H.plot(X_LV_L, H_LV_L, color=lc_LV_L, lw=lw)
        H_LV_V=H_PTX(P0+TT*0, T_LV_V, X_LV_V)
        ax_H.plot(X_LV_V, H_LV_V, color=lc_LV_V, lw=lw)
        ax_H.fill(np.append(np.append(np.flip(X_LV_V),X_LV_L),Xmin), 
                np.append(np.append(np.flip(H_LV_V),H_LV_L),prop1.H_v), 
                fc=fc_LV, ec='None', lw=lw, alpha=alpha, label='L+V')
        # log scale branch
        ax_H_log.fill(np.append(X_LV_L, np.flip(X_log)), np.append(H_LV_L, np.flip(H_log_VLH1)),fc=fc_LV, ec='None', lw=lw, alpha=alpha, label='L+V')
        ax_H_log.plot(X_LV_L, H_LV_L, color=lc_LV_L, lw=lw)
        ax_H_log.plot(X_LV_V, H_LV_V, color=lc_LV_V, lw=lw)
        ax_H_log.fill_betweenx(H_LV_V, X_LV_V, Xmin,fc='w', ec='None', lw=lw, alpha=alpha, zorder=zorder_p-1)
        # print(H_LV_L,H_LV_V)
        # second part
        TT=T_LV_L[1]
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(TT,TT*0+P0)))
        ax_T.plot(X_LV_L, TT, color=lc_LV_L)
        ax_T.fill_betweenx(TT, X_LV_L, Xmin, fc=fc_LV, ec='None', alpha=alpha)
        H_LV_L=H_PTX(P0+TT*0, TT, X_LV_L)
        ax_H.plot(X_LV_L, H_LV_L, color=lc_LV_L)
        ax_H.fill_between(np.insert(X_LV_L,0,Xmin), np.insert(H_LV_L,0,prop2.H_v), Hmax, fc=fc_LV, ec='None', alpha=alpha)
        # log scale branch
        ax_H_log.fill_between(X_log, H_log_VLH2, Hmax, fc=fc_LV, ec='None', alpha=alpha)

    else:
        TT=T_LH[0]
        X_LH=np.array(sw.Mol2Wt(sw.X_HaliteLiquidus(TT, TT*0+P0)))
        ax_T.plot(X_LH, TT, color=lc_LH)
        ax_T.fill_betweenx(TT, X_LH, Tmax, fc=fc_LH, ec='None', alpha=alpha,label='L+H')
        H_LH=H_PTX(P0+TT*0, TT, X_LH)
        ax_H.plot(X_LH, H_LH, color=lc_LH)
        H_LH_bottom=H_PTX(P0+X_LH*0, Tmin+X_LH*0, X_LH)
        ax_H.fill_between(X_LH, H_LH, H_LH_bottom,fc=fc_LH, ec='None', alpha=alpha,label='L+H')
        # LV
        TT=T_LV_L[0]
        X_LV_L=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_LiquidBranch(TT,TT*0+P0)))
        X_LV_L[0]=sw.Mol2Wt(X_crit_P0)
        ax_T.plot(X_LV_L, TT, color=lc_LV_L)
        ax_T.fill_between(X_LV_L, TT, Tmax, fc=fc_LV, ec='None', alpha=alpha, label='L+V')
        X_LV_V=np.array(sw.Mol2Wt(sw.X_VaporLiquidCoexistSurface_VaporBranch(T_LV_V,T_LV_V*0+P0)))
        X_LV_V[0]=sw.Mol2Wt(X_crit_P0)
        ax_T.plot(X_LV_V, T_LV_V, color=lc_LV_V)
        ax_T.fill_between(X_LV_V, T_LV_V, Tmax, fc=fc_LV, ec='None', alpha=alpha)
        H_LV_L=H_PTX(P0+TT*0, TT, X_LV_L)
        ax_H.plot(X_LV_L, H_LV_L, color=lc_LV_L)
        ax_H.fill_between(X_LV_L, H_LV_L, Hmax, fc=fc_LV, ec='None', alpha=alpha, label='L+V')
        H_LV_V=H_PTX(P0+TT*0, T_LV_V, X_LV_V)
        ax_H.plot(X_LV_V, H_LV_V, color=lc_LV_V)
        ax_H.fill_between(X_LV_V, H_LV_V, Hmax, fc=fc_LV, ec='None', alpha=alpha)
    # critical point
    ax_T.scatter(sw.Mol2Wt(X_crit_P0), T_crit_P0, fc='k',ec='r', zorder=zorder_p, label='Critical point')
    ax_H.scatter(sw.Mol2Wt(X_crit_P0), prop_crit.H, fc='k',ec='r', zorder=zorder_p, label='Critical point')
    if(sw.Mol2Wt(X_crit_P0)<=X_max_log):
        ax_H_log.scatter(sw.Mol2Wt(X_crit_P0), prop_crit.H, fc='k',ec='r', label='Critical point',zorder=zorder_p)
    # print(sw.Mol2Wt(X_crit_P0), prop_crit.H)
    # print(sw.Mol2Wt(X_crit_P0), T_crit_P0)
    # valid region
    fill_invalid=ax_H.fill_between(X, H_boundary_Tmax, Hmax, fc='gray', ec='k', label='T>%.0f $^{\circ}$C'%(Tmax))
    ax_H.fill_between(X, H_boundary_Tmin, 0, fc='lightgray', ec='k',label='T<%.0f $^{\circ}$C'%(Tmin))
    ax_H_log.fill_between(X, H_boundary_Tmin, 0, fc='lightgray', ec='k',label='T<=%.0f $^{\circ}$C'%(Tmin))
    ax_H.axhspan(1,2,fc='w',ec='k',label='L')
    ax_T.axhspan(1,2,fc='w',ec='k',label='L')
    ax_T.legend(ncol=2)
    ax_H.legend(ncol=2)
    # plot phase state with H0, and bulk composition
    prop0=sw.prop_pHX(P0*1E5, H0, X0)
    if(savefig==True):
        print('H=%.3f MJ/kg, X=%.1f wt%%NaCl, P=%.1f bar:\n  T=%.1f C\n  H_l=%.3f MJ/kg, H_v=%.3f MJ/kg, H_h=%.3f MJ/kg\n  X_l=%.1f, X_v=%.2E\n  S_l=%.1f, S_v=%.1f, S_h=%.1f'%(H0/1E6, X0*100, P0, prop0.T, prop0.H_l/1E6, prop0.H_v/1E6, prop0.H_h/1E6, prop0.X_l*100, prop0.X_v*100, prop0.S_l*100, prop0.S_v*100, prop0.S_h*100))
    # plot in X-H space
    ms_base=15
    ms_l=ms_base*prop0.S_l
    if(ms_l<ms_base/4):
        ms_l=ms_base/4
    ms_v=10*prop0.S_v
    if(ms_v<ms_base/4):
        ms_v=ms_base/4
    if(prop0.X_l<=X_max_log):
        ax_H_log.plot(prop0.X_l, prop0.H_l,'o',mfc='r',mec='c',markersize=ms_l, zorder=10,alpha=1)
    else:
        ax_H.plot(prop0.X_l, prop0.H_l,'o',mfc='r',mec='c',markersize=ms_l, zorder=10,alpha=1)
    if(prop0.X_v>=X_max_log):
        ax_H.plot(prop0.X_v, prop0.H_v,'o',mfc='b',mec='y', markersize=ms_v, zorder=10,alpha=1)
    else:
        ax_H_log.plot(prop0.X_v, prop0.H_v,'o',mfc='b',mec='y', markersize=ms_v, zorder=10,alpha=1)
    # plot in X-T space
    ax_T.plot(prop0.X_l, prop0.T,'o',mfc='r',mec='c',markersize=ms_l, zorder=10,alpha=1)
    if(prop0.X_v>=Xmin):
        ax_T.plot(prop0.X_v, prop0.T,'o',mfc='b',mec='y', markersize=ms_v, zorder=10,clip_on=False,alpha=1)
    ax_T.text(0.02,0.1, 'T=%.0f $^{\circ}$C'%(prop0.T),ha='left',va='bottom',color='k',fontweight='bold', transform=ax_T.transAxes)
    # axis
    ax_H.xaxis.set_major_locator(MultipleLocator(0.2))
    ax_H.xaxis.set_minor_locator(MultipleLocator(0.04))
    for ax in axes:
        if(ax==ax_T):
            ax.set_xlim(Xmin, Xmax)
        else:
            ax.set_xlim(X_max_log, Xmax)
        ax.set_xlabel('Salinity (wt% NaCl)')
        ticks=ax.get_xticks()
        ax.xaxis.set_ticks(ticks)
        labels=[]
        for tick in ticks:
            labels.append(str('%.0f'%(tick*100)))
        ax.xaxis.set_ticklabels(labels)
    ax_H.set_xlim(X_max_log, 1)
    ax_T.set_ylim(Tmin,Tmax)
    ax_T.set_ylabel('Temperature ($^{\circ}$C)')
    ax_H.set_ylabel('Specific enthalpy (MJ/kg)')
    # enthalpy axis
    ticks=ax_H.get_yticks()
    ax_H.yaxis.set_ticks(ticks)
    labels=[]
    for tick in ticks:
        labels.append(str('%.1f'%(tick/1E6)))
    ax_H.yaxis.set_ticklabels(labels)
    ax_H.set_ylim(0,Hmax)
    ax_H_log.set_ylim(ax_H.get_ylim())
    # write pressure value
    ax_H_log.text(0.01,0.01,'P=%.0f MPa'%(P0/10),ha='left',va='bottom',transform=ax_H_log.transAxes,zorder=10)
    ax_T.text(0.01,0.01,'P=%.0f MPa'%(P0/10),ha='left',va='bottom',transform=ax_T.transAxes,zorder=10)

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.98)
    figname='PhaseDiagram_XT_P0_%.0fbar'%(P0)
    if(savefig):
        for fmt in fmt_figs:
            plt.savefig('%s.%s'%(figname,fmt), dpi=dpi)
    else:
        return ax_T, ax_H, ax_H_log, figname, prop0
    # plt.show()
# Heat a fluid with initial composition(bulk salinity) by increasing bulk enthalpy,
# visualize the phase change path
def makeAnimation_heatingSaltWater(P0=250,fmt='png',dpi=400,figpath='.',X_max_log=0.1):
    H=np.arange(0.1,4.0,0.01)*1E6
    H_l, H_v, X_l, X_v,T_l,T_v=[],[],[],[],[],[]

    for i in range(0,len(H)):
        H0=H[i]
        ax_T, ax_H, ax_H_log, figname,prop=XT_P0_logLinearScale(P0=P0,X0=0.032,H0=H0,savefig=False,X_max_log=X_max_log)
        if(prop.X_l>=1E-6):
            X_l.append(prop.X_l)
            H_l.append(prop.H_l)
            T_l.append(prop.T)
        if(prop.X_v>=1E-6):
            X_v.append(prop.X_v)
            H_v.append(prop.H_v)
            T_v.append(prop.T)
        # plot history path
        # X-H
        ax_H_log.plot(X_l, H_l, '.', color='r',zorder=11,alpha=1,ms=3)
        ax_H.plot(X_l, H_l, '.', color='r',zorder=11, clip_on=False,alpha=1,ms=3)
        
        ax_H.plot(X_v, H_v, '.', color='b',zorder=11, clip_on=False,alpha=1,ms=3)
        ax_H_log.plot(X_v, H_v, '.', color='b',zorder=11,alpha=1,ms=3)
        # X-T
        ax_T.plot(X_v,T_v,'.',color='b',zorder=11, clip_on=False,alpha=1,ms=3)
        ax_T.plot(X_l,T_l,'.',color='r',zorder=11, clip_on=False,alpha=1,ms=3)
        plt.savefig('%s/%04d.%s'%(figpath,i,fmt),dpi=dpi)
        print('Progress: %d/%d'%(i, len(H)))
    # XT_P0_logLinearScale(P0=P0,H0=2.65*1E6)
def main(argv):
    # argc=len(argv)
    # usage(argv)
    # exit(0)
    XT_P0()
    # makeAnimation_heatingSaltWater(figpath='animation')
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))