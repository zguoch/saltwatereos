#!/Users/zguo/.pyenv/shims/python
# -*-coding:utf-8-*-
# Calculate specific enthalpy of halite
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/04/06, GEOMAR
# ===============================================================
import sys
import os
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
print(pyswEOS.__version__)
from pyswEOS import H2ONaCl
from pyswEOS import H2O
from pyswEOS import NaCl
water=H2O.cH2O()
sw=H2ONaCl.cH2ONaCl()
halite=NaCl.cNaCl()

def writeInputFile_sowt(T_c, p_bar, x_wt,filename):
    fpout=open(filename,'w')
    for i in range(0,len(T_c)):
        fpout.write('%f\n%f\n%f\n\n'%(T_c[i],p_bar[i],sw.Wt2Mol(x_wt[i])))
    fpout.write('0\n100\n0\n\n')
    fpout.close()
def SpecificEnthalpy_halite(T_c,p_bar):
    H=np.zeros_like(T_c)
    for i in range(0,len(T_c)):
        H[i]=halite.SpecificEnthalpy(T_c[i],p_bar[i])
    return H
def Cp_halite(T_c,p_bar):
    Cp=np.zeros_like(T_c)
    for i in range(0,len(T_c)):
        Cp[i]=halite.Cp(T_c[i],p_bar[i])
    return Cp
def SpecificEnthalpy_sw(T_c,p_bar,x_wt):
    H=np.zeros_like(T_c)
    for i in range(0,len(T_c)):
        prop=sw.prop_pTX(p_bar[i]*1E5, T_c[i]+273.15, x_wt[i])
        H[i]=prop.H
    return H
def cal_Plot():
    P_max_LVH  = H2ONaCl.P_max_LVH
    T_Pmax_LVH = H2ONaCl.T_Pmax_LVH
    T     = np.linspace(T_Pmax_LVH-0.1, T_Pmax_LVH+0.1, 10)
    P     = P_max_LVH + np.linspace(-2,2,3)
    # plot
    fig=plt.figure(figsize=(10,6))
    ax=plt.gca()
    # ax_cp=ax.twinx()
    for p0 in P:
        H=SpecificEnthalpy_halite(T, p0+T*0)/1000  # kJ/kg
        # Cp=Cp_halite(T, p0+T*0)
        H_sw=SpecificEnthalpy_sw(T, p0+T*0, 1+T*0)/1000
        # H
        l,=ax.plot(T,H,marker='.',label='%.8f bar'%(p0))
        ax.plot(T,H_sw,'o',mfc='None',mec='r')
        # # Cp
        # ax_cp.plot(T,Cp,color=l.get_color(),ls='-.',marker='^')
        # write input file for sowat_ptx
        writeInputFile_sowt(T,p0+T*0, 1+T*0,'P%.4fbar_halite.txt'%(p0))
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('Specific enthalpy (kJ/kg)')
    ax.grid(which='major',axis='both',lw=0.1,color='k')
    ax.grid(which='minor',axis='both',lw=0.05,color='gray')
    ax.legend()
    ax.set_title('Specific enthalpy of halite close to LVH region')
    plt.tight_layout()
    for fmt in ['pdf','svg','jpg']:
        plt.savefig('Halite_H_Cp.%s'%(fmt))
    # plt.show()
def main(argv):
    cal_Plot()
if __name__ == '__main__':
    sys.exit(main(sys.argv))