#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Plot benchmark test result of swEOS
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2020/11/15, GEOMAR
# ===============================================================

import sys
import argparse
import sys
import os
from colored import fg, bg, attr
C_GREEN = fg('green')
C_RED = fg('red')
C_BLUE = fg('blue')
C_DEFAULT = attr('reset')
#===============================================================
import linecache
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math

def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Plot benchmark test result of swEOS'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2020/11/15, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))
datapath='build'
figpath='../doxygen/images'

mH2O = 18.015268
mNaCl = 58.4428
T_Triple_NaCl = 800.7
P_Triple_NaCl = 0.0005
w_singleFig=6
def XmolPercent2XwtPercent(molPercent):
    Tmp = molPercent / 100.0
    wtPercent = mNaCl * Tmp / (mNaCl * Tmp + (1 - Tmp) * mH2O) * 100
    return wtPercent
def plot_HaliteLiquidus(fname0='X_HaliteLiquidus',fmt='svg'):
    fig,axs=plt.subplots(1,2,figsize=(w_singleFig*2,5))
    ax=axs[0]
    # Fig 7a
    P=[500, 2000, 4000] # bar
    linestyles=['dotted','dashed','solid']
    for p,ls in zip(P,linestyles):
        fname=str('%s/%s_P%.0fbar.dat'%(datapath,fname0,p))
        data=np.loadtxt(fname)
        T=data[:,0]
        X=data[:,1]
        ax.plot(X,T,ls=ls,label='%.0f bar'%(p))
    ax.set_xlim(0,1)
    ax.set_ylim(0,1000)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(0.04))
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.grid(which='major',color='gray',lw=0.03)
    ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('X$_{\mathregular{NaCl}}$ [Mole fraction of NaCl]')
    ax.set_ylabel('T [$^{\circ}$C]')
    ax.legend(loc='lower right')
    # wt% NaCl
    ax2=ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.xaxis.set_major_locator(ax.xaxis.get_major_locator())
    ax2.xaxis.set_minor_locator(ax.xaxis.get_minor_locator())
    ticklabels_wt=[]
    for tick in ax2.get_xticks():
        ticklabels_wt.append('%.1f'%(XmolPercent2XwtPercent(tick*100)))
    ax2.xaxis.set_ticklabels(ticklabels_wt)
    ax2.set_xlabel('X$_{\mathregular{NaCl}}$ [wt. % NaCl]')
    ax.text(0.02,0.98,'(a)',transform=ax.transAxes,va='top',ha='left',fontweight='bold')
    # Fig 7b
    ax=axs[1]
    arryT=[25] # degC
    linestyles=['solid','dashed','dotted']
    for T,ls in zip(arryT,linestyles):
        fname=str('%s/%s_T%.0fC.dat'%(datapath,fname0,T))
        data=np.loadtxt(fname)
        P=data[:,0]
        X=data[:,1]
        ax.plot(X,P,ls=ls,label='%.0f $^{\circ}$C'%(T))
    ax.set_xlim(0.1,0.115)
    ax.set_ylim(0,5000)
    ax.xaxis.set_major_locator(MultipleLocator(0.005))
    ax.yaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(0.001))
    ax.yaxis.set_minor_locator(MultipleLocator(200))
    ax.grid(which='major',color='gray',lw=0.03)
    ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('X$_{\mathregular{NaCl}}$ [Mole fraction of NaCl]')
    ax.set_ylabel('P [bar]')
    ax.legend(loc='lower right')
    ax.text(0.02,0.98,'(b)',transform=ax.transAxes,va='top',ha='left',fontweight='bold')
    # wt% NaCl
    ax2=ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.xaxis.set_major_locator(ax.xaxis.get_major_locator())
    ax2.xaxis.set_minor_locator(ax.xaxis.get_minor_locator())
    ticklabels_wt=[]
    for tick in ax2.get_xticks():
        ticklabels_wt.append('%.1f'%(XmolPercent2XwtPercent(tick*100)))
    ax2.xaxis.set_ticklabels(ticklabels_wt)
    ax2.set_xlabel('X$_{\mathregular{NaCl}}$ [wt. % NaCl]')

    figname=str('%s/%s.%s'%(figpath,'HaliteLiquidus',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_HaliteMelting(fname0='HaliteMeltingCurve.dat',fmt='svg'):
    data=np.loadtxt('%s/%s'%(datapath,fname0))
    T=data[:,0]
    P=data[:,1]
    fig=plt.figure(figsize=(w_singleFig,w_singleFig))
    ax=plt.gca()
    ax.plot(T,P)
    ax.set_xlim(800,930)
    ax.set_ylim(0,5000)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(200))
    ax.grid(which='major',color='gray',lw=0.03)
    ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('T [$^{\circ}$C]')
    ax.set_ylabel('P [bar]')

    figname=str('%s/%s.%s'%(figpath,'HaliteMeltingCurve',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_SublimationBoiling(fname0='HaliteSublimationBoilingCurve.dat',fmt='svg'):
    data=np.loadtxt('%s/%s'%(datapath,fname0))
    T=data[:,0]
    P_subl=data[:,1]
    P_boil=data[:,2]
    fig=plt.figure(figsize=(w_singleFig,w_singleFig))
    ax=plt.gca()
    ax.semilogy(T,P_subl,label='Sublimation curve')
    ax.semilogy(T,P_boil,label='Boiling curve')
    # triple point
    ax.semilogy(T_Triple_NaCl,P_Triple_NaCl,'o',mfc='red',mec='w',label='Triple point of NaCl')
    ax.semilogy([T_Triple_NaCl, T_Triple_NaCl],[ax.get_ylim()[1],P_Triple_NaCl],ls='dashed')
    # text of phase region
    ax.text(T_Triple_NaCl-50,10**-1.5,'Halite',va='top',ha='right',fontweight='bold')
    ax.text(T_Triple_NaCl+50,10**-1.5,'Liquid',va='top',ha='left',fontweight='bold')
    ax.text(T_Triple_NaCl,1E-5,'Vapor',va='top',ha='left',fontweight='bold')
    ax.text(T_Triple_NaCl+10,P_Triple_NaCl,'(%.2f, %.2E)'%(T_Triple_NaCl,P_Triple_NaCl),va='top',ha='left',color='r')
    ax.set_xlim(300, 1100)
    ax.set_ylim(1E-14, 1E-1)
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=20))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8),numticks=20))
    ax.grid(which='major',color='gray',lw=0.03)
    ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('T [$^{\circ}$C]')
    ax.set_ylabel('P [bar]')
    ax.legend(loc='lower right')
    

    figname=str('%s/%s.%s'%(figpath,'HaliteSublimationBoilingCurves',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_CriticalPressure_Salinity(fname0='HaliteCritical_P_X',fmt='svg'):
    fig,axs=plt.subplots(2,2,figsize=(w_singleFig*2,8),sharex='col',gridspec_kw={'wspace':0.05,'hspace':0.05})
    # full range 
    data=np.loadtxt('%s/%s_fullrange.dat'%(datapath,fname0))
    T=data[:,0]
    P=data[:,1]
    X=data[:,2]
    # Pressure
    ax=axs[0][0]
    ax.plot(T,P)
    ax.set_ylim(0,2500)
    ax.set_ylabel('P [bar]')
    ax.yaxis.set_major_locator(MultipleLocator(500))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
    ax.text(0.02,0.98,'(a)',va='top',ha='left',fontweight='bold',transform=ax.transAxes)
    # Salinity
    ax=axs[1][0]
    ax.plot(T,X)
    ax.set_ylim(0,0.15)
    ax.set_ylabel('X$_{\mathregular{NaCl}}$ [[Mole fraction of NaCl]')
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.text(0.02,0.98,'(c)',va='top',ha='left',fontweight='bold',transform=ax.transAxes)
    ax.set_xlim(T.min(),T.max())
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))

    # close to water critical temperature
    data=np.loadtxt('%s/%s_closeWaterCriticalT.dat'%(datapath,fname0))
    T=data[:,0]
    P=data[:,1]
    X=data[:,2]
    # Pressure
    ax=axs[0][1]
    ax.plot(T,P)
    ax.set_ylim(200, 310)
    ax.set_ylabel('P [bar]')
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(2))
    ax.text(0.02,0.98,'(d)',va='top',ha='left',fontweight='bold',transform=ax.transAxes)
    # Salinity
    ax=axs[1][1]
    ax.plot(T,X)
    ax.set_ylim(0,0.03)
    ax.set_ylabel('X$_{\mathregular{NaCl}}$ [[Mole fraction of NaCl]')
    ax.yaxis.set_major_locator(MultipleLocator(0.005))
    ax.yaxis.set_minor_locator(MultipleLocator(0.001))
    ax.text(0.02,0.98,'(c)',va='top',ha='left',fontweight='bold',transform=ax.transAxes)
    ax.set_xlim(T.min(),T.max())
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))

    for i in range(0,2):
        axs[i][1].yaxis.set_ticks_position('right')
        axs[i][1].yaxis.set_label_position('right')
        axs[1][i].set_xlabel('T [$^{\circ}$C]')
        for j in range(0,2):
            axs[i][j].grid(which='major',color='gray',lw=0.03)
            axs[i][j].grid(which='minor',color='lightgray',lw=0.03)

    figname=str('%s/%s.%s'%(figpath,'HaliteCriticalCurves',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_HaliteSaturatedVaporComposition(fname0='X_HaliteSaturatedVapor',fmt='svg'):
    fig,axs=plt.subplots(1,2,figsize=(w_singleFig*2,5))
    ax=axs[1]
    # Fig 8b
    P=[1, 10, 50, 100, 200,300] # bar
    for p in P:
        fname=str('%s/%s_P%.0fbar.dat'%(datapath,fname0,p))
        data=np.loadtxt(fname)
        T=data[:,0]
        X=data[:,1]
        ax.semilogx(X,T,label='%.0f bar'%(p))
    # plot X_VLH
    fname=str('%s/X_VLH.dat'%(datapath))
    data=np.loadtxt(fname)
    T=data[:,0]
    X=data[:,1]
    ax.semilogx(X,T,label='V+L+H')
    
    ax.set_xlim(1E-12, 1)
    ax.set_ylim(100, 850)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=20))
    ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.1, 0.2,0.3, 0.4,0.5, 0.6,0.7, 0.8, 0.9),numticks=20))
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.grid(which='major',color='gray',lw=0.03)
    ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('log$_{\mathregular{10}}$(X$_{\mathregular{NaCl}}$) [Mole fraction of NaCl]')
    ax.set_ylabel('T [$^{\circ}$C]')
    ax.legend(loc='lower right',ncol=2)
    ax.text(0.02,0.98,'(b)',transform=ax.transAxes,va='top',ha='left',fontweight='bold')
    # Fig 8a
    ax=axs[0]
    arryT=[450, 500, 550] # degC
    linestyles=['solid','dashed','dotted']
    for T in arryT:
        fname=str('%s/%s_T%.0fC.dat'%(datapath,fname0,T))
        data=np.loadtxt(fname)
        P=data[:,0]
        X=data[:,1]
        l,=ax.semilogx(X,P,label='%.0f $^{\circ}$C'%(T))
        l_Fig8,=ax.semilogx(X/10,P,ls='dotted',color=l.get_color())
    ax.set_xlim(1E-12, 1E-4)
    ax.set_ylim(0,400)
    # ax.xaxis.set_major_locator(MultipleLocator(0.005))
    ax.yaxis.set_major_locator(MultipleLocator(50))
    # ax.xaxis.set_minor_locator(MultipleLocator(0.001))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.grid(which='major',color='gray',lw=0.03)
    ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('log$_{\mathregular{10}}$(X$_{\mathregular{NaCl}}$) [Mole fraction of NaCl]')
    ax.set_ylabel('P [bar]')
    leg1=ax.legend(loc='lower right')
    ax.add_artist(leg1)
    ax.legend(handles=(l,l_Fig8),loc='center left',labels=['Original','Fig. 8 of Driesner & Heinrich(2007)\nX$_{\mathregular{NaCl}}$ is shifted 10 times to the left'])
    ax.vlines(2E-5,ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1],ls='dashed',color='gray')
    ax.text(0.02,0.98,'(a)',transform=ax.transAxes,va='top',ha='left',fontweight='bold')
    
    figname=str('%s/%s.%s'%(figpath,'HaliteSaturatedVaporComposition',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_P_VLH(fname0='P_VLH',fmt='svg'):
    fig,axs=plt.subplots(1,2,figsize=(w_singleFig*2,5),gridspec_kw={'wspace':0.05})
    # full range 
    data=np.loadtxt('%s/%s_fullrange.dat'%(datapath,fname0))
    T=data[:,0]
    P=data[:,1]
    # Pressure
    ax=axs[0]
    ax.plot(T,P)
    ax.set_ylim(0,400)
    ax.set_ylabel('P [bar]')
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.text(0.02,0.98,'(a)',va='top',ha='left',fontweight='bold',transform=ax.transAxes)
    ax.set_xlim(T.min(),T.max())
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))

    # low temperature
    data=np.loadtxt('%s/%s_lowT.dat'%(datapath,fname0))
    T=data[:,0]
    P=data[:,1]
    ax=axs[1]
    ax.semilogy(T,P)
    ax.set_ylim(1E-2, 1E2)
    ax.set_ylabel('log$_{\mathregular{10}}$(P [bar])')
    ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=20))
    ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.1, 0.2,0.3, 0.4,0.5, 0.6,0.7, 0.8, 0.9),numticks=20))
    ax.text(0.02,0.98,'(b)',va='top',ha='left',fontweight='bold',transform=ax.transAxes)
    ax.set_xlim(T.min(),T.max())
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))

    axs[1].yaxis.set_ticks_position('right')
    axs[1].yaxis.set_label_position('right')
    for i in range(0,2):
        axs[i].set_xlabel('T [$^{\circ}$C]')
        axs[i].grid(which='major',color='gray',lw=0.03)
        axs[i].grid(which='minor',color='lightgray',lw=0.03)

    figname=str('%s/%s.%s'%(figpath,'Pressure_VLH',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_X_VL(fname0='X_VL',fmt='svg'):
    rows,cols=4,3
    fig,axs=plt.subplots(rows,cols,figsize=(w_singleFig*cols,4*rows),gridspec_kw={'wspace':0.15,'hspace':0.2})
    # Fig 12
    arryT=np.array([200, 300, 350, 375, 380, 400, 500, 600, 800, 1000]) # degC, 375 is 375.5
    label_subfig=['a','b','c','d','e','f','g','h','i','j','k','l']
    xmin_axes=np.array([-12, -8, -7, -5, -5, -5, -4, -4, -5, -4])
    major_locator_y=[1, 10, 20, 20, 20, 50, 100, 100, 500, 500]
    linestyles=['solid','dashed','dotted']
    for i in range(0,rows):
        for j in range(0,cols):
            ax=axs[i][j]
            index = j+i*cols
            if(index>(len(arryT)-1)):
                ax.axis('off')
                continue
            T=arryT[index]
            fname=str('%s/%s_LiquidBranch_T%.0fC.dat'%(datapath,fname0,T))
            fname_vaporBranch=str('%s/%s_VaporBranch_T%.0fC.dat'%(datapath,fname0,T))
            data=np.loadtxt(fname)
            data_vaporBranch=np.loadtxt(fname_vaporBranch)
            P=data[:,0]
            X=data[:,1]
            P_vaporBranch=data_vaporBranch[:,0]
            X_vaporBranch=data_vaporBranch[:,1]
            l_l,=ax.semilogx(X,P,label='Liquid branch')
            l_v,=ax.semilogx(X_vaporBranch,P_vaporBranch,label='Vapor branch')
            ax.set_xlim(10.0**xmin_axes[index], 1)
            ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=20))
            ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.1, 0.2,0.3, 0.4,0.5, 0.6,0.7, 0.8, 0.9),numticks=20))
            ax.yaxis.set_major_locator(MultipleLocator(major_locator_y[index]))
            ax.yaxis.set_minor_locator(MultipleLocator(major_locator_y[index]/5.0))
            ax.grid(which='major',color='gray',lw=0.03)
            ax.grid(which='minor',color='lightgray',lw=0.03)
            ax.set_xlabel('log$_{\mathregular{10}}$(X$_{\mathregular{NaCl}}$) [Mole fraction of NaCl]')
            ax.set_ylabel('P [bar]')
            ax.text(0.02,0.98,'(%s)'%label_subfig[index],transform=ax.transAxes,va='top',ha='left',fontweight='bold')
            str_T='%.0f $^{\circ}$C'%arryT[index]
            if(arryT[index]==375):
                str_T = '375.5 $^{\circ}$C'
            ax.text(0.5,0.2,str_T,transform=ax.transAxes,va='bottom',ha='center', fontsize=12,fontweight='bold')
            ax.legend(loc='lower left')
    figname=str('%s/%s.%s'%(figpath,'X_VaporLiquidCoexistSurface',fmt))
    plt.savefig(figname, bbox_inches='tight')
def plot_water_boilingCurve(fname0='Water_P_boiling.dat',fmt='svg'):
    data=np.loadtxt('%s/%s'%(datapath,fname0))
    T=data[:,0]
    P1=data[:,1]
    P2=data[:,2]
    rho_v_sat=data[:,3]
    rho_l_sat=data[:,4]
    fig=plt.figure(figsize=(w_singleFig,w_singleFig))
    ax=plt.gca()
    ax.plot(T[0:-1],P1[0:-1],label='cH2O::P_Boiling')
    ax.plot(T[0:-1],P2[0:-1],label='cH2O::BoilingCurve')
    ax.plot(T[-1],P1[-1],'o',mfc='red',mec='w',label='Critical point')
    ax2=ax.twinx()
    ax2.plot(T[0:-1],rho_v_sat[0:-1])
    ax2.plot(T[0:-1],rho_l_sat[0:-1])
    ax.legend()
    # ax.set_xlim(800,930)
    # ax.set_ylim(0,5000)
    # ax.xaxis.set_major_locator(MultipleLocator(50))
    # ax.yaxis.set_major_locator(MultipleLocator(1000))
    # ax.xaxis.set_minor_locator(MultipleLocator(10))
    # ax.yaxis.set_minor_locator(MultipleLocator(200))
    # ax.grid(which='major',color='gray',lw=0.03)
    # ax.grid(which='minor',color='lightgray',lw=0.03)
    ax.set_xlabel('T [$^{\circ}$C]')
    ax.set_ylabel('P [bar]')

    figname=str('%s/%s.%s'%(figpath,'water_boilingCurve',fmt))
    plt.savefig(figname, bbox_inches='tight')
def main(argv):
    # argc=len(argv)
    # usage(argv)
    # exit(0)

    # plot_HaliteMelting()
    # plot_SublimationBoiling()
    # plot_CriticalPressure_Salinity()
    # plot_HaliteLiquidus()
    # plot_HaliteSaturatedVaporComposition()
    # plot_P_VLH()
    # plot_X_VL()
    plot_water_boilingCurve()

if __name__ == '__main__':
    sys.exit(main(sys.argv))