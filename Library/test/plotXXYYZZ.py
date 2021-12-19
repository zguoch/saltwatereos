#!python
# -*-coding:utf-8-*-
# plot result of lookup
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/12/19, GEOMAR
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
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math

def usage(argv):
    basename = argv[0].split('/')
    basename = basename[len(basename)-1]
    description='plot result of lookup'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/12/19, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' xx.txt yy.txt zz.txt'+C_DEFAULT)
    print('='*len(head))


def main(argv):
    if(len(argv)!=4):
        usage(argv)
        exit()
    
    xx = np.loadtxt(argv[1])
    yy = np.loadtxt(argv[2])
    zz = np.loadtxt(argv[3])
    fig = plt.figure(figsize=(8,8))
    ax=plt.gca()
    ax.contourf(xx,yy,zz, levels=50, cmap='rainbow')
    plt.savefig('result.pdf', bbox_inches='tight')


if __name__ == '__main__':
    sys.exit(main(sys.argv))