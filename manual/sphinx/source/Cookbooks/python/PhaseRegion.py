import pyswEOS
from pyswEOS import H2ONaCl
sw=H2ONaCl.cH2ONaCl()

p=200 #bar
T=400 #deg.C
X=0.032 #wt.% NaCl
region = sw.findPhaseRegion(T,p,X)
region_name=sw.getPhaseRegionName(region);
print(" Pressure(bar): ",p)
print(" Temperature(deg.C): ",T)
print(" Salinity (wt.% NaCl): ",X)
print(" Phase region index: ",region)
print(" Phase region name: ",region_name)