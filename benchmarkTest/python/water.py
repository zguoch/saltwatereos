#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Test water (H2O) sub functions
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
    description='Test water (H2O) sub functions'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/05/24, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))

def T_P_Boiling():
    fig=plt.figure()
    ax=plt.gca()

    T=np.linspace(H2O.TMIN,H2O.T_Critic,100)
    P=np.zeros_like(T)
    for i in range(0,len(T)):
        P[i]=water.P_Boiling(T[i])
    ax.plot(T,P,label='Given T')
    P=np.linspace(1, H2O.P_Critic, 50)
    T=np.zeros_like(P)
    for i in range(0,len(P)):
        T[i]=water.T_Boiling(P[i])
    ax.plot(T,P,'.',label='Given P')
    ax.set_xlabel('T ($^{\circ}$C)')
    ax.set_ylabel('P (bar)')
    ax.legend()
    ax.set_title('Boiling curve of water')
    plt.show()
def main(argv):
    T_P_Boiling()

if __name__ == '__main__':
    sys.exit(main(sys.argv))