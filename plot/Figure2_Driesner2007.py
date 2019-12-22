import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd

def plot_figure2(datapath, filename, islog=False):
    fname=datapath+filename+'.csv'
    data=pd.read_csv(fname)
    xlabel=''
    if(islog):
        X=data['X'].values
        xlabel='log$_{10}$(wt. % NaCl)'
    else:
        X=data['X'].values*100
        xlabel='wt. % NaCl'
    P=data['P'].values/1e5
    T=data['T'].values
    reg=data['Region'].values
    reg_unique=np.unique(reg)
    # calculate mesh
    nP=0
    nX=0
    if(P[0]==P[1]):
        nP=1
        for i in range(1,len(P)):
            if(P[i]==P[i-1]):
                nP=nP+1
            else:
                break
        nX=int(len(P)/nP)
    else:
        nX=1
        for i in range(1,len(X)):
            if(X[i]==X[i-1]):
                nX=nX+1
            else:
                break
        nP=int(len(X)/nX)
    shape0=(nP,nX)
    P=np.reshape(P,shape0)
    X=np.reshape(X,shape0)
    reg=np.reshape(reg,shape0)
    reg_unique=np.unique(reg)
    print(reg_unique)

    plt.figure()
    ax=plt.gca()

    # scatter=ax.scatter(X,P,c=reg,s=0.1, alpha=0.5)
    cf=ax.contourf(X,P,reg,cmap='Dark2', vmin=0,vmax=7)
    ax.set_xlim(np.min(X),np.max(X))
    ax.set_ylim(np.min(P),np.max(P))
    ax.set_ylabel('P[bar]')
    ax.set_xlabel(xlabel)
    ax.text(np.min(X)+0.03*(np.max(X)-np.min(X)),np.max(P)-0.1*(np.max(P)-np.min(P)),str(T[0])+' $^{\circ}$C',fontsize=12)
    # plt.colorbar(cf)
    plt.savefig('c++/'+filename+'.pdf')
    # plt.show()

datapath='../build/'

# figure 2a
plot_figure2(datapath,'Figure2a_Driesner2007')
# figure 2c
plot_figure2(datapath,'Figure2c_Driesner2007')
# figure 2e
plot_figure2(datapath,'Figure2e_Driesner2007')
# figure 2g
plot_figure2(datapath,'Figure2g_Driesner2007')
# figure 2i
plot_figure2(datapath,'Figure2i_Driesner2007')
# figure 2k
plot_figure2(datapath,'Figure2k_Driesner2007')

# log xaxis
# figure 2b
plot_figure2(datapath,'Figure2b_Driesner2007',islog=True)
# figure 2d
plot_figure2(datapath,'Figure2d_Driesner2007',islog=True)
# figure 2f
plot_figure2(datapath,'Figure2f_Driesner2007',islog=True)
# figure 2h
plot_figure2(datapath,'Figure2h_Driesner2007',islog=True)
# figure 2j
plot_figure2(datapath,'Figure2j_Driesner2007',islog=True)