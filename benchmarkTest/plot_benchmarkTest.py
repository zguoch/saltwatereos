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
from matplotlib.ticker import MultipleLocator
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
def XmolPercent2XwtPercent(molPercent):
    Tmp = molPercent / 100.0
    wtPercent = mNaCl * Tmp / (mNaCl * Tmp + (1 - Tmp) * mH2O) * 100
    return wtPercent
def plot_HaliteLiquidus(fname0='X_HaliteLiquidus',fmt='svg'):
    fig,axs=plt.subplots(1,2,figsize=(12,5))
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
    ax.set_ylabel('T[$^{\circ}$C]')
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

def main(argv):
    # argc=len(argv)
    # usage(argv)
    # exit(0)
    plot_HaliteLiquidus()

if __name__ == '__main__':
    sys.exit(main(sys.argv))