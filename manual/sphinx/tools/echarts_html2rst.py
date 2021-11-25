#!python
# -*-coding:utf-8-*-
# Extract elements in body tag of html generated by pyecharts to rst file
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'Zhikui Guo, 2021/11/19, GEOMAR
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
    description='Extract elements in body tag of html generated by pyecharts to rst file'
    num_symbol=int((len(description)+20 - len(basename))/2)
    head='='*num_symbol+basename+'='*num_symbol
    print(head)
    print(description)
    print('Zhikui Guo, 2021/11/19, GEOMAR')
    print('[Example]: '+C_RED + basename+C_BLUE + ' example usage'+C_DEFAULT)
    print('='*len(head))

def echarts_html2rst(file_html):
    file_rst=file_html.replace('.html','.rst')
    fpout=open(file_rst,'w')
    alldata=linecache.getlines(file_html)
    start,end=0,0
    for i in range(0,len(alldata)):
        line=alldata[i]
        if(line[0:6]=='<body>'):
            start = i 
        elif(line[0:7]=='</body>'):
            end = i
    fpout.write('.. raw:: html\n\n')
    for i in range(start, end+1):
        # line = alldata[i].replace('\n','')
        fpout.write('\t%s'%(alldata[i]))
    fpout.close()
def echarts_html2rst_dir(htmlDir):
    for file in os.listdir(htmlDir):
        if(file.split('.')[-1]=='html'):
            echarts_html2rst('%s/%s'%(htmlDir,file))
            print('%s processed'%(file.replace('.html','')))
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--save', help='save animation',action='store_true')
    parser.add_argument('input',type=str, help='Input grd file')
    parser.add_argument('-o','--output',type=str, help='output path',default='.')
    args = parser.parse_args()
    input=args.input
    if(os.path.isdir(input)):
        echarts_html2rst_dir(input)
    elif(os.path.isfile(input)):
        if(input.split('.')[-1]=='html'):
            echarts_html2rst(input)
        else:
            print('File is not html: %s\n'%(input))


if __name__ == '__main__':
    sys.exit(main(sys.argv))