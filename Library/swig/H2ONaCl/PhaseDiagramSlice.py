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
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
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
fmt_figs=['pdf','svg']
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
# which='x', 'y', value=xvalue, yvalue
def getManualClabel(CS,levels_c, which, value):
    manual_locations=[]
    # specify clabel position
    if(which=='x'):
        x_clabel=value # approximate x of clabel
        manual_locations=np.zeros((len(levels_c),2))
        for i in range(0,len(CS.collections)):
            dist_min=np.zeros((len(CS.collections[i]._paths),1))
            for j in range(0,len(CS.collections[i]._paths)):
                _path=CS.collections[i]._paths[j]
                vertices_cs=_path._vertices
                dist_x=np.abs(vertices_cs[:,0]-x_clabel)
                dist_min[j]=dist_x.min()
            ind_path=np.argmin(dist_min)
            vertices_cs=CS.collections[i]._paths[ind_path]._vertices
            dist_x=np.abs(vertices_cs[:,0]-x_clabel)
            ind=np.argmin(dist_x)
            manual_locations[i,:]=vertices_cs[ind]
    elif(which=='y'):
        y_clabel=value # approximate y of clabel
        manual_locations=np.zeros((len(levels_c),2))
        for i in range(0,len(CS.collections)):
            dist_min=np.zeros((len(CS.collections[i]._paths),1))
            for j in range(0,len(CS.collections[i]._paths)):
                _path=CS.collections[i]._paths[j]
                vertices_cs=_path._vertices
                dist_y=np.abs(vertices_cs[:,1]-y_clabel)
                dist_min[j]=dist_y.min()
            ind_path=np.argmin(dist_min)
            vertices_cs=CS.collections[i]._paths[ind_path]._vertices
            dist_y=np.abs(vertices_cs[:,1]-y_clabel)
            ind=np.argmin(dist_y)
            manual_locations[i,:]=vertices_cs[ind]
    else:
        print('only support which="x" or "y"')
    return manual_locations
def H_PTX(P_bar, T_C, X_wt):
    H=np.zeros_like(P_bar)
    for j in range(0,len(H)):
        prop=sw.prop_pTX(P_bar[j]*1E5, T_C[j]+273.15, X_wt[j])
        H[j]=prop.H
    return H
def XT_P0(P0=300, Xmin=1E-5, Xmax=0.999999, Tmin=1, Tmax=1000, nX=100, nT=100, 
    fc_VLH=(1,1,0,0.8), ec_VLH='b',fc_VH=(0.5,0,0.5,0.5), lc_LH='b',lc_LV_L='b',
    lc_LV_V='r',fc_LH='lightblue',fc_LV='lightgreen',datapath='.'):
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
    T_nearCrit=np.arange(T_crit_P0, T_crit_P0+2, 0.01)
    T_LH, T_LV_V = [np.linspace(Tmin, Tmax, 100)], np.append(T_nearCrit, np.linspace(T_crit_P0+2, Tmax, 100))
    T_LV_L=[T_LV_V]
    # plot
    # fig=plt.figure()
    fig,axes=plt.subplots(1,2,gridspec_kw={"width_ratios":[1,1],'wspace':0.05},figsize=(14,6))
    ax_T,ax_H=axes[0],axes[1]
    ax_H.yaxis.set_ticks_position('right')
    ax_H.yaxis.set_label_position('right')
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
    ax_H.set_ylim(0,Hmax)
    ticks=ax_H.get_yticks()
    ax_H.yaxis.set_ticks(ticks)
    labels=[]
    for tick in ticks:
        labels.append(str('%.1f'%(tick/1E6)))
    ax_H.yaxis.set_ticklabels(labels)
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
        x_tmp,y_tmp=np.array([Xmin, HX[2], Xmax, Xmin]),np.array([prop1.H_v, prop1.H_l, prop1.H_h, prop1.H_v])
        ax_H.fill(x_tmp,y_tmp, fc=fc_VLH,ec=ec_VLH, label='L+V+H')
        xy_display_V, xy_display_H = ax_H.transData.transform((0,prop1.H_v)), ax_H.transData.transform((1,prop1.H_h))
        dxdy_display=np.abs(xy_display_V-xy_display_H)
        angle_V_H=np.arctan2(dxdy_display[1],dxdy_display[0])/np.pi*180
        ax_H.text(HX[2]*1.05, prop1.H_l*1.05,'%.1f$^{\circ}$C'%(prop1.T),fontsize=9,va='center',ha='center',rotation=-angle_V_H)
        ax_H.fill([Xmin, HX[3], Xmax, Xmin],[prop2.H_v, prop2.H_l, prop2.H_h, prop2.H_v], fc=fc_VLH,ec=ec_VLH)
        xy_display_V, xy_display_H = ax_H.transData.transform((0,prop2.H_v)), ax_H.transData.transform((1,prop2.H_h))
        dxdy_display=np.abs(xy_display_V-xy_display_H)
        angle_V_H=np.arctan2(dxdy_display[1],dxdy_display[0])/np.pi*180
        ax_H.text(HX[3]*0.97, prop2.H_l*0.97,'%.1f$^{\circ}$C'%(prop2.T),fontsize=9,va='center',ha='center',rotation=-angle_V_H)
        
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
        T_LV_L=[np.linspace(T_crit_P0, prop1.T-1, len(T_LV_V)), np.linspace(prop2.T+1, Tmax, len(T_LV_V))]
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
    fill_invalid=ax_H.fill_between(X, H_boundary_Tmax, Hmax, fc='gray', ec='k', label='T>%.0f $^{\circ}$C'%(Tmax))
    ax_H.fill_between(X, H_boundary_Tmin, 0, fc='lightgray', ec='k',label='T<%.0f $^{\circ}$C'%(Tmin))
    ax_H.axhspan(1,2,fc='w',ec='b',label='L')
    ax_T.axhspan(1,2,fc='w',ec='b',label='L')
    ax_T.legend(ncol=2)
    ax_H.legend(ncol=2)
    
    # ax_H.grid(axis='both',which='both')

    # calculate temperature in X-H-P space
    # X=np.linspace(1E-5,1,200)
    H=np.linspace(1E-5,4.5,200)*1E6
    XX_XH,HH_XH=np.meshgrid(X,H)
    TT_XH=np.zeros_like(XX_XH)
    T=np.linspace(Tmin,Tmax,400)
    XX_XT,TT_XT=np.meshgrid(X,T)
    HH_XT=np.zeros_like(XX_XT)
    Rho_XT=np.zeros_like(XX_XT)
    dataname_T='T_p%.0fMPa_H_X'%(P0/10)
    dataname_H='H_p%.0fMPa_H_X'%(P0/10)
    dataname_Rho='Rho_p%.0fMPa_H_X'%(P0/10)
    datafile_T='%s/%s.dat'%(datapath,dataname_T)
    datafile_H='%s/%s.dat'%(datapath,dataname_H)
    datafile_Rho = '%s/%s.dat'%(datapath,dataname_Rho)
    calculate=False
    fpout_T,fpout_H=None,None
    if(os.path.exists(datafile_T)):
        data=np.loadtxt(datafile_T)
        if(data.shape==TT_XH.shape):
            TT_XH=np.loadtxt(datafile_T)
            HH_XT=np.loadtxt(datafile_H)
            Rho_XT=np.loadtxt(datafile_Rho)
        else:
            calculate=True
            fpout_T=open(datafile_T,'w')
            fpout_H=open(datafile_H,'w')
            fpout_Rho=open(datafile_Rho,'w')
    else:
        calculate=True
        fpout_T=open(datafile_T,'w')
        fpout_H=open(datafile_H,'w')
        fpout_Rho=open(datafile_Rho,'w')
    if(calculate):
        for i in range(0,len(H)):
            for j in range(0,len(X)):
                if(H[i]>H_boundary_Tmax[j]):
                    TT_XH[i][j]=np.nan
                else:
                    prop=sw.prop_pHX(P0*1E5, H[i], X[j])
                    TT_XH[i][j]=prop.T
                fpout_T.write('%.6E '%(TT_XH[i][j]))
            fpout_T.write('\n')
            print('Progress %d/%d'%(i,len(H)))
        fpout_T.close()
        for i in range(0,XX_XT.shape[0]):
            for j in range(0,XX_XT.shape[1]):
                prop=sw.prop_pTX(P0*1E5, TT_XT[i,j]+273.15, XX_XT[i,j])
                fpout_H.write('%.6E '%(prop.H))
                fpout_Rho.write('%.6E '%(prop.Rho))
            fpout_H.write('\n')
            fpout_Rho.write('\n')
            print('Progress %d/%d'%(i,len(H)))
        fpout_H.close()
        fpout_Rho.close()
    # plot Temperature contour
    levels_c=np.array([500, 600, 700, 800, 900])
    CS=ax_H.contour(XX_XH,HH_XH,TT_XH,levels=levels_c,colors='k',linewidths=0.5)
    ax_H.clabel(CS,fmt='%.0f $^{\circ}$C',fontsize=9,manual=getManualClabel(CS,levels_c,'x',0.2))
    levels_c=np.array([3,100,200,300,350])
    CS=ax_H.contour(XX_XH,HH_XH,TT_XH,levels=levels_c,colors='k',linewidths=0.5)
    ax_H.clabel(CS, inline=True, fmt='%.0f $^{\circ}$C', fontsize=9,manual=getManualClabel(CS,levels_c,'x',0.23))

    levels_c=np.array([390, 400, 420, 440])
    CS=ax_H.contour(XX_XH,HH_XH,TT_XH,levels=levels_c,colors='k',linewidths=0.5, linestyles='dashed')
    ax_H.clabel(CS, inline=True, fmt='%.0f $^{\circ}$C', fontsize=9,manual=getManualClabel(CS,levels_c,'y',1.95E6))
    # plot H in X-T space
    levels_c=np.arange(1.0,3.4,0.2)
    CS=ax_T.contour(XX_XT, TT_XT, HH_XT/1E6, levels=levels_c,colors=['r','k','r','k','r','k','r','k','r','k'],linewidths=0.5)
    ax_H.clabel(CS, inline=True, fmt='%.1f MJ/kg', fontsize=9,manual=getManualClabel(CS,levels_c,'y',600))
    ax_H.clabel(CS, inline=True, fmt='%.1f MJ/kg', fontsize=9,manual=getManualClabel(CS,levels_c,'y',800))
    levels_c=np.arange(3.4,4,0.2)
    CS=ax_T.contour(XX_XT, TT_XT, HH_XT/1E6, levels=levels_c,colors=['r','k','r','k','r','k','r','k','r','k'],linewidths=0.5)
    ax_H.clabel(CS, inline=True, fmt='%.1f MJ/kg', fontsize=9,manual=getManualClabel(CS,levels_c,'y',800))
    levels_c=np.arange(0.2,1.0,0.2)
    CS=ax_T.contour(XX_XT, TT_XT, HH_XT/1E6, levels=levels_c,colors=['r','k','r','k','r','k','r','k','r','k'],linewidths=0.5)
    ax_H.clabel(CS, inline=True, fmt='%.1f MJ/kg', fontsize=9,manual=getManualClabel(CS,levels_c,'x',0.5))
    # ---plot Rho
    ax_T.contourf(XX_XT, TT_XT, Rho_XT, levels=50, cmap='rainbow')
    # ---------------------------------
    

    # write pressure value
    ax_H.text(0.01,0.02,'P=%.0f MPa'%(P0/10),ha='left',va='bottom',transform=ax_H.transAxes,zorder=10)
    ax_T.text(0.01,0.02,'P=%.0f MPa'%(P0/10),ha='left',va='bottom',transform=ax_T.transAxes,zorder=10)

    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.98)
    for fmt in fmt_figs:
        plt.savefig('PhaseDiagram_XT_P0_%.0fbar.%s'%(P0,fmt), dpi=dpi)
    # plt.show()
line_intersection = lambda p1, p2, x0: (x0-p1[0])/(p2[0]-p1[0])*(p2[1]-p1[1])+p1[1]
def XT_P0_logLinearScale(P0=150, X0=None, H0=None, Xmin=0.5E-4, Xmax=0.999999, Tmin=5, Tmax=1000, nX=100, nT=100, 
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
    T_nearCrit=np.arange(T_crit_P0, T_crit_P0+2, 0.01)
    T_LH, T_LV_V = [np.linspace(Tmin, Tmax, 100)], np.append(T_nearCrit, np.linspace(T_crit_P0+2, Tmax, 100))
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
    ax_H_log.spines['right'].set_visible(False)
    ax_H.spines['left'].set_visible(False)
    ax_H_log.set_yticks([])
    ax_H_log.axvspan(0,1,fc='lightgray',alpha=0.4,zorder=10)
    # ax_H_log.text(0.5,0.01, 'Log scale',ha='center',va='bottom',color='b', transform=ax_H_log.transAxes)
    ax_H_log.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=20))
    ax_H_log.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8),numticks=20))
    ax_H_log.tick_params(axis='x',which='both', colors='blue')
    ax_H_log.spines['bottom'].set_color('b')
    # enthalpy axis
    ticks=ax_H.get_yticks()
    ax_H.yaxis.set_ticks(ticks)
    labels=[]
    for tick in ticks:
        labels.append(str('%.1f'%(tick/1E6)))
    ax_H.yaxis.set_ticklabels(labels)
    ax_H.set_ylim(0,Hmax)
    ax_H_log.set_ylim(ax_H.get_ylim())
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
    ax_H_log.set_xlim((Xmin, X_max_log))
    ax_T.set_ylim(Tmin,Tmax)
    ax_T.set_ylabel('Temperature ($^{\circ}$C)')
    ax_H.set_ylabel('Specific enthalpy (MJ/kg)')
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
        T_LV_L=[np.linspace(T_crit_P0, prop1.T-1, len(T_LV_V)), np.linspace(prop2.T+1, Tmax, len(T_LV_V))]
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
        ind_valid_VL_V=(H_LV_V<=H_log_VLH1.max())
        x_tmp=np.append(np.append(np.flip(X_LV_V[ind_valid_VL_V]), X_LV_L),np.flip(X_log))
        y_tmp=np.append(np.append(np.flip(H_LV_V[ind_valid_VL_V]), H_LV_L),np.flip(H_log_VLH1))
        ax_H_log.fill(x_tmp,y_tmp,fc=fc_LV, ec='None', lw=lw, alpha=alpha, label='L+V')
        # ax_H_log.fill(np.append(X_LV_L, np.flip(X_log)), np.append(H_LV_L, np.flip(H_log_VLH1)),fc=fc_LV, ec='None', lw=lw, alpha=alpha, label='L+V')
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
        # log scale
        ax_H_log.fill_between(X_LV_V, H_LV_V, Hmax, fc=fc_LV, ec='None', alpha=alpha)
        ax_H_log.plot(X_LV_V, H_LV_V, color=lc_LV_V)
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
    if(not ((H0==None) | (X0==None))):
        ax_H_log.text(0.01,0.06, 'H=%.2f MJ/kg'%(H0/1E6),ha='left',va='bottom',color='k',fontweight='bold', transform=ax_H_log.transAxes,zorder=10)
        prop0=sw.prop_pHX(P0*1E5, H0, X0)
        rhom=prop0.S_l*prop0.Rho_l + prop0.S_v*prop0.Rho_v + prop0.S_h*prop0.Rho_h
        X_f = (prop0.S_l*prop0.Rho_l*prop0.X_l + prop0.S_v*prop0.Rho_v*prop0.X_v)/rhom
        X_s = (prop0.S_h*prop0.Rho_h)/rhom
        Xm = X_f + X_s
        if(savefig==True):
            print('H=%.3f MJ/kg, X=%.1f wt%%NaCl, P=%.1f bar:\n  T=%.1f C\n  H_l=%.3f MJ/kg, H_v=%.3f MJ/kg, H_h=%.3f MJ/kg\n  X_l=%.1f, X_v=%.2E\n  S_l=%.1f, S_v=%.1f, S_h=%.1f'%(H0/1E6, X0*100, P0, prop0.T, prop0.H_l/1E6, prop0.H_v/1E6, prop0.H_h/1E6, prop0.X_l*100, prop0.X_v*100, prop0.S_l*100, prop0.S_v*100, prop0.S_h*100))
            Xm=(prop0.S_l*prop0.Rho_l*prop0.X_l + prop0.S_v*prop0.Rho_v*prop0.X_v + prop0.S_h*prop0.Rho_h)/rhom
            print('  Bulk density: %.2f\n  Bulk salinity: %.2f\n  Fluid salinity: %.2f'%(rhom, X_f + X_s, X_f))
        # plot in X-H space
        c_l,c_v,c_h,c_f='r','b','k','c'
        ms_base=15
        ms_l=ms_base*prop0.S_l
        ms_v=ms_base*prop0.S_v
        ms_h=ms_base*prop0.S_h
        if(ms_l<ms_base/4):
            ms_l=ms_base/4
        if(ms_v<ms_base/4):
            ms_v=ms_base/4
        if((prop0.S_h>1E-6) & (ms_h<ms_base/4)):
            ms_h=ms_base/4
        ec_l,ec_v,ec_h,fc_l,fc_v,fc_h=c_v, c_l,c_h,c_l,c_v,c_h
        if((prop0.S_l<1E-6) & (prop0.S_h>1E-6)):  #V+H
            ec_v,ec_h=c_h,c_v
        if((prop0.S_h>1E-6) & (prop0.S_v>1E-6) & (prop0.S_l>1E-6)): #V+L+H
            ms_h=ms_base*(X_s/Xm)
        if(prop0.X_l<=X_max_log):
            ax_H_log.plot(prop0.X_l, prop0.H_l,'o',mfc=fc_l,mec=ec_l,markersize=ms_l, zorder=10,alpha=1)
        else:
            ax_H.plot(prop0.X_l, prop0.H_l,'o',mfc=fc_l,mec=ec_l,markersize=ms_l, zorder=10,alpha=1)
        if(prop0.X_v>=X_max_log):
            ax_H.plot(prop0.X_v, prop0.H_v,'o',mfc=fc_v,mec=ec_v, markersize=ms_v, zorder=10,alpha=1)
        else:
            ax_H_log.plot(prop0.X_v, prop0.H_v,'o',mfc=fc_v,mec=ec_v, markersize=ms_v, zorder=10,alpha=1)
        ax_H.plot(1, prop0.H_h,'o',mfc=fc_h,mec=ec_h, markersize=ms_h, zorder=10,alpha=1, clip_on=False) # Halite
        # plot bulk composition in multiphase region
        if((prop0.S_h>1E-6) | (prop0.S_v>1E-6)):
            ax_bulk_base=ax_H
            if(X0<=X_max_log):
                ax_bulk_base=ax_H_log
            min_display,max_display,data_display=ax_bulk_base.transAxes.transform((0,0)),ax_bulk_base.transAxes.transform((1,1)),ax_bulk_base.transData.transform((X0,H0))
            x0y0_frac=(data_display-min_display)/(max_display-min_display)
            w_ax,h_ax=1, 0.07
            ax_bulkH_log = ax_bulk_base.inset_axes([x0y0_frac[0]-w_ax/2, x0y0_frac[1]-h_ax/2, w_ax,h_ax], transform=ax_bulk_base.transAxes, zorder=10)
            r1,r2=1,0.4
            ax_bulkH_log.pie([prop0.S_l, prop0.S_v, prop0.S_h],colors=['r','b','k'], shadow=False,radius=r1, wedgeprops=dict(width=r2, edgecolor='w',lw=0.1, alpha=0.8),startangle=0)
            ax_bulkH_log.pie([X_f/Xm, X_s/Xm],colors=['c','k'], shadow=False,radius=r1-r2, wedgeprops=dict(width=r1-r2-0.05, edgecolor='w',lw=0.1, alpha=0.8),startangle=0)
            ax_bulkH_log.xaxis.set_ticks([])
            ax_bulkH_log.yaxis.set_ticks([])
            ax_bulkH_log.axis('scaled')
            # ax_bulk_base.plot(X0, H0, 'o', mfc='gray',mec='None',markersize=ms_base/2,zorder=11,alpha=0.5)
            
        # plot phase state and saturation
        w_rect,h_rect,x0,y0=0.8,0.04,0.01,0.12
        ax_H_log.text(x0, y0+h_rect*1.1,'Phase saturation (%)\n',ha='left',va='bottom',color='gray',fontweight='normal',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        if(prop0.S_l>1E-6):
            ax_H_log.text(x0, y0+h_rect*1.1,'Liquid',ha='left',va='bottom',color=c_l,fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        if(prop0.S_v>1E-6):
            ax_H_log.text(x0+0.3, y0+h_rect*1.1,'Vapor',ha='left',va='bottom',color=c_v,fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        if(prop0.S_h>1E-6):
            ax_H_log.text(x0+0.6, y0+h_rect*1.1,'Halite',ha='left',va='bottom',color=c_h,fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        # phase saturation
        x0_tmp,w_tmp=x0,w_rect*prop0.S_l
        rect_l = mpatches.Rectangle((x0_tmp, y0),w_tmp,h_rect,fc='r', ec="None",transform=ax_H_log.transAxes,zorder=zorder_p)
        ax_H_log.add_patch(rect_l)
        if(prop0.S_l>0.15):
            ax_H_log.text(x0_tmp+w_tmp/2, y0+h_rect/2,'%.0f%%'%(prop0.S_l*100),ha='center',va='center',color='w',fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        x0_tmp,w_tmp=x0_tmp+w_tmp,w_rect*prop0.S_v
        rect_v = mpatches.Rectangle((x0_tmp, y0),w_tmp,h_rect,fc='b', ec="None",transform=ax_H_log.transAxes,zorder=zorder_p)
        ax_H_log.add_patch(rect_v)
        if(prop0.S_v>0.15):
            ax_H_log.text(x0_tmp+w_tmp/2, y0+h_rect/2,'%.0f%%'%(prop0.S_v*100),ha='center',va='center',color='w',fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        x0_tmp,w_tmp=x0_tmp+w_tmp,w_rect*prop0.S_h
        rect_h = mpatches.Rectangle((x0_tmp, y0),w_tmp,h_rect,fc='k', ec="None",transform=ax_H_log.transAxes,zorder=zorder_p)
        ax_H_log.add_patch(rect_h)
        if(prop0.S_h>0.15):
            ax_H_log.text(x0_tmp+w_tmp/2, y0+h_rect/2,'%.0f%%'%(prop0.S_h*100),ha='center',va='center',color='w',fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        # salinity
        w_rect,h_rect,x0,y0=0.8,0.05,0.01,0.25
        ax_H_log.text(x0, y0+h_rect*1.1,'Bulk composition (wt.% NaCl)\n',ha='left',va='bottom',color='gray',fontweight='normal',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        ax_H_log.text(x0, y0+h_rect*1.1,'Fluid',ha='left',va='bottom',color=c_f,fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        if(prop0.S_h>1E-6):
            ax_H_log.text(x0+w_rect, y0+h_rect*1.1,'Solid (100% NaCl)',ha='right',va='bottom',color=c_h,fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        
        x0_tmp,w_tmp=x0,w_rect*(X_f/Xm)
        rect_f = mpatches.Rectangle((x0_tmp, y0),w_tmp,h_rect,fc=c_f, ec="None",transform=ax_H_log.transAxes,zorder=zorder_p)
        ax_H_log.add_patch(rect_f)
        if((X_f/Xm)>0.15):
            ax_H_log.text(x0_tmp+w_tmp/2, y0+h_rect/2,'%.0f%%'%(X_f*100),ha='center',va='center',color='w',fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        x0_tmp,w_tmp=x0_tmp+w_tmp,w_rect*(X_s/Xm)
        rect_s = mpatches.Rectangle((x0_tmp, y0),w_tmp,h_rect,fc=c_h, ec="None",transform=ax_H_log.transAxes,zorder=zorder_p)
        ax_H_log.add_patch(rect_s)
        # if((X_s/Xm)>0.15):
        #     ax_H_log.text(x0_tmp+w_tmp/2, y0+h_rect/2,'%.0f%%'%(X_s*100),ha='center',va='center',color='w',fontweight='bold',fontsize=10, transform=ax_H_log.transAxes,zorder=10)
        
        # x0_tmp,w_tmp=x0_tmp+w_tmp,w_rect*prop0.S_h
        # rect_h = mpatches.Rectangle((x0_tmp, y0),w_tmp,h_rect,fc='k', ec="None",transform=ax_H_log.transAxes)
        # ax_H_log.add_patch(rect_h)

        # plot in X-T space
        ax_T.plot(prop0.X_l, prop0.T,'o',mfc=fc_l,mec=ec_l,markersize=ms_l, zorder=10,alpha=1)
        if(prop0.X_v>=Xmin):
            ax_T.plot(prop0.X_v, prop0.T,'o',mfc=fc_v,mec=ec_v, markersize=ms_v, zorder=10,clip_on=False,alpha=1)
        ax_T.text(0.02,0.1, 'T=%.0f $^{\circ}$C'%(prop0.T),ha='left',va='bottom',color='k',fontweight='bold', transform=ax_T.transAxes)
    
    # write pressure value
    ax_H_log.text(0.01,0.02,'P=%.0f MPa'%(P0/10),ha='left',va='bottom',transform=ax_H_log.transAxes,zorder=10)
    ax_T.text(0.01,0.02,'P=%.0f MPa'%(P0/10),ha='left',va='bottom',transform=ax_T.transAxes,zorder=10)

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
def makeAnimation_heatingSaltWater(P0=250,X0=0.2,fmt='png',dpi=400,figpath='.',X_max_log=0.1):
    # calculate VLH points
    HminHmaxXminXmax=sw.HX_VaporLiquidHaliteCoexist(P0) 
    HX=HminHmaxXminXmax
    H=[]
    if(len(HX)==4):
        prop_VLH1=sw.prop_pHX(P0*1E5, HX[0], HX[2]) # first VLH point
        prop_VLH2=sw.prop_pHX(P0*1E5, HX[1], HX[3]) #second VLH point
        prop1_LVH1=sw.prop_pTX(P0*1E5, prop_VLH1.T+273.15-1, X0)
        prop2_LVH1=sw.prop_pTX(P0*1E5, prop_VLH1.T+273.15+1, X0)
        prop1_LVH2=sw.prop_pTX(P0*1E5, prop_VLH2.T+273.15-1, X0)
        prop2_LVH2=sw.prop_pTX(P0*1E5, prop_VLH2.T+273.15+1, X0)
        prop1_LVH1.H,prop2_LVH1.H=prop1_LVH1.H/1E6,prop2_LVH1.H/1E6
        prop1_LVH2.H,prop2_LVH2.H=prop1_LVH2.H/1E6,prop2_LVH2.H/1E6
        H=np.arange(1.3,prop1_LVH1.H,0.02)*1E6  # L -> L+V
        H_lvh=np.arange(prop1_LVH1.H,prop2_LVH1.H,0.001)*1E6 # LVH, refine
        H=np.append(H, H_lvh)
        H_vh=np.arange(prop2_LVH1.H,prop1_LVH2.H,0.01)*1E6 # VH
        H=np.append(H, H_vh)
        H_lvh=np.arange(prop1_LVH2.H,prop2_LVH2.H,0.001)*1E6 # LVH, refine
        H=np.append(H, H_lvh)
        H_vl=np.arange(prop2_LVH2.H,3.40,0.01)*1E6 # LV
        H=np.append(H, H_vl)
    else:
        H=np.arange(0.1,3.8,0.02)*1E6
    
    H_l, H_v, X_l, X_v,T_l,T_v=[],[],[],[],[],[]

    for i in range(0,len(H)):
        H0=H[i]
        ax_T, ax_H, ax_H_log, figname,prop=XT_P0_logLinearScale(P0=P0,X0=X0,H0=H0,savefig=False,X_max_log=X_max_log)
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
        os.system('cp %s/%04d.%s %s/%04d.%s'%(figpath,i,fmt, figpath,len(H)*2-1-i,fmt))
        print('Progress: %d/%d'%(i, len(H)))
    # XT_P0_logLinearScale(P0=P0,H0=1.6*1E6, X0=0.2)
def main(argv):
    # argc=len(argv)
    # usage(argv)
    # exit(0)
    # XT_P0(P0=10)
    # XT_P0(P0=15)
    # XT_P0(P0=150)
    XT_P0(P0=250)
    # XT_P0(P0=350)
    # XT_P0(P0=450)
    # XT_P0(P0=200)
    # XT_P0(P0=300)
    # XT_P0(P0=500)
    # makeAnimation_heatingSaltWater(figpath='animation')
    # XT_P0_logLinearScale(P0=650,H0=1.6*1E6, X0=0.2)
    # XT_P0_logLinearScale(P0=250,H0=1.6*1E6, X0=0.2)
    # XT_P0_logLinearScale(P0=650)

    

if __name__ == '__main__':
    sys.exit(main(sys.argv))