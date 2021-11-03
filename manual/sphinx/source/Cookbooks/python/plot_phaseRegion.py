import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import pyswEOS
from pyswEOS import H2ONaCl
from pyswEOS import H2O
from pyswEOS import NaCl
water=H2O.cH2O()
sw=H2ONaCl.cH2ONaCl()
halite=NaCl.cNaCl()
def getPhaseRegion(TT,PP,XX):
    region=np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            print(sw.findPhaseRegion(TT[i][j],PP[i][j],XX[i][j])[0])
            region[i][j] = sw.findPhaseRegion(TT[i][j],PP[i][j],XX[i][j])[0]
    return region
def plotRegion(ax,x,y,region,showColorbar=False,xlabel='',ylabel='',title='',label=''):
    minlevel,maxlevel=0,9
    nlevel=int(maxlevel-minlevel)+1
    norm = mpl.colors.BoundaryNorm(np.linspace(minlevel-0.5,maxlevel+0.5, nlevel+1), nlevel) # 划重点
    CSf=ax.contourf(x,y,region+0.5,cmap=plt.get_cmap("Paired", nlevel), norm=norm,vmin=0,vmax=7)
    if(showColorbar):
        cbar_kw=dict(ticks=np.arange(minlevel,maxlevel+1))
        ax_cb=ax.inset_axes([1.02,0,0.03,1],transform=ax.transAxes)
        cbar = plt.colorbar(CSf, cax=ax_cb,label='Phase region',**cbar_kw)
        cbar.set_ticks(np.array([0,1,2,3,4,5,6,7])+0.5)
        cbar.set_ticklabels(['L','VL','V','L+H','V+H','VLH','VL-L','VL-V'])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.text(0.02,0.98,label,va='top',ha='left',transform=ax.transAxes,fontsize=14,color='w')
def plotRegions(Tmin=H2ONaCl.TMIN_C,Tmax=H2ONaCl.TMAX_C, nT=300,
                Pmin=H2ONaCl.PMIN, Pmax=H2ONaCl.PMAX, nP=300, 
                Xmin=H2ONaCl.XMIN, Xmax=H2ONaCl.XMAX, nX=300):
    # figure
    w_singleFig=6
    fig,axs=plt.subplots(1,3,figsize=(w_singleFig*3,4),gridspec_kw={'wspace':0.2,'hspace':0.05})
    X_label,T_label,P_label='Bulk salinity (wt% NaCl)','Temperature ($^{\circ}$C)','Pressure (bar)'
    # constant P
    P0=400 #bar
    X=np.linspace(Xmin,Xmax,nX)
    T=np.linspace(Tmin,Tmax,nT)
    TT,XX=np.meshgrid(T,X)
    PP=P0+TT*0
    region_P0=getPhaseRegion(TT,PP,XX)
    ax=axs[0]
    plotRegion(ax,XX*100,TT,region_P0,xlabel=X_label,ylabel=T_label,title='Pressure=%.0f bar'%(P0),label='(a)')

    # constant T
    T0=400 #deg.C
    X=np.linspace(Xmin,Xmax,nX)
    P=np.linspace(Pmin,800,nP)
    PP,XX=np.meshgrid(P,X)
    TT=T0+PP*0
    region_T0=getPhaseRegion(TT,PP,XX)
    ax=axs[1]
    plotRegion(ax,XX,PP,region_T0,xlabel=X_label,ylabel=P_label,title='Temperature=%.0f $^{\circ}$C'%(T0),label='(b)')

    # constant X
    X0=0.032 #deg.C
    T=np.linspace(Tmin,Tmax,nT)
    P=np.linspace(Pmin,Pmax,nP)
    TT,PP=np.meshgrid(T,P)
    XX=X0+PP*0
    region_X0=getPhaseRegion(TT,PP,XX)
    ax=axs[2]
    plotRegion(ax,TT,PP,region_X0,xlabel=T_label,ylabel=P_label,title='Salinity=%.2f wt%% NaCl'%(X0*100),showColorbar=True,label='(c)')

    plt.show()
plotRegions()
 