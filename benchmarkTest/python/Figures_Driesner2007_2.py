#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Fig. 5 of Driesner(2007) Part 2.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/05/24, GEOMAR
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
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
import linecache
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import string
import pyswEOS
print('pyswEOS version: ',pyswEOS.__version__)
from pyswEOS import H2ONaCl
from pyswEOS import H2O
from pyswEOS import NaCl
water=H2O.cH2O()
sw=H2ONaCl.cH2ONaCl()
halite=NaCl.cNaCl()
figpath='.'
fmt_figs = ['jpg','svg','pdf']
def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='Fig. 5 of Driesner(2007) Part 2.'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/05/24, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))
def molV_PX(P, X):
    """Calculate molar volume by giving pressure and salinity.

    Args:
        P (float): pressure [bar]
        X (float): salinity [wt.% NaCl]
    """
    arrayT=np.linspace(1,1000,1000)
    V_water=np.zeros_like(arrayT)
    X_mol=sw.Wt2Mol(X/100)
    for i in range(0,len(arrayT)):
        T=arrayT[i]
        T_star = sw.T_star_V(T, P, X_mol) # deg.C
        V_water[i] = H2O.MolarMass / water.Rho(T_star, P) # m^3/mol
    return {'T':arrayT, 'V_water':V_water}
def Fig2():
    P0 = 1000 #bar
    T0 = 600 # deg.C
    fig = plt.figure()
    ax = plt.gca()
    ax.set_title('p=%.0f bar'%(P0))
    for X in [0,10]:
        data = molV_PX(P0,X)
        ax.plot(data['T'],data['V_water']*1E6,label='%.0f wt.%% NaCl'%(X))
    ax.set_ylabel('V ($cm^3/mol$)')
    ax.set_xlabel('T ($^{\circ}$C)')
    ax.set_xlim(200,700)
    ax.set_ylim(15,50)
    ax.legend()
    # plot illustration of T-T^star
    T_star = sw.T_star_V(T0, P0, sw.Wt2Mol(10/100)) # deg.C
    V_star = H2O.MolarMass / water.Rho(T_star, P0)*1E6
    ax.plot(T_star,V_star,'o',mfc='w',mec='b',zorder=10)
    ax.plot(T0,V_star,'o',mfc='w',mec='r',zorder=10)
    ax.axhline(V_star,lw=0.5,ls='dotted',color='k')
    l, = ax.plot([T0,T0],[ax.get_ylim()[0], V_star],ls='dashed',color='r',lw=0.5)
    l_star, = ax.plot([T_star,T_star],[ax.get_ylim()[0], V_star],ls='dashed',color='b',lw=0.5)
    ax.text(T_star*0.98,ax.get_ylim()[0],'T$^{*}$',ha='right',va='bottom',color=l_star.get_color())
    ax.text(T0*1.02,ax.get_ylim()[0],'T',ha='left',va='bottom',color=l.get_color())
    for fmt in fmt_figs:
        figname=str('%s/Fig2_Driesner2007_part2.%s'%(figpath,fmt))
        plt.savefig(figname, bbox_inches='tight')
def main(argv):
    # argc=len(argv)
    # usage(argv)
    # exit(0)
    Fig2()

if __name__ == '__main__':
    sys.exit(main(sys.argv))