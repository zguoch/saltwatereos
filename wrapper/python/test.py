import linecache
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
import H2O

water = H2O.IAPWS95()

# T=np.linspace(0,400,100)
# P=np.zeros_like(T)
# for i in range(0,len(T)):
#     P[i]=water.P_Boiling(T[i])
T=np.linspace(0, 1000, 200)
P=np.linspace(1, 1000, 200)
PP,TT=np.meshgrid(P,T)
rho=np.zeros_like(TT)
for i in range(0,len(T)):
    for j in range(0,len(P)):
        rho[i][j]=water.Rho(T[i],P[j])
fig=plt.figure()
ax=plt.gca()
# ax.plot(T,P)
ax.contourf(TT,PP,rho,nlevel=50)
plt.show()