import H2ONaCl
import H2O
water=H2O.cH2O()
sw=H2ONaCl.cH2ONaCl()
# 1. function
print(sw.Rho_brine(300, 200, 0.2))
print(sw.P_VaporLiquidHaliteCoexist(200))
# 2. function with multi returns
p_crit,x_crit=sw.P_X_Critical(800)
print(p_crit,x_crit,sw.Mol2Wt(x_crit)*100)
# 3. constants
print(H2ONaCl.TMAX)
# 4. 
print(water.P_Boiling(200))

# 5. calculate properties
prop=sw.prop_pHX(100E5, 2E6, 0.5)
print(prop.T,sw.getPhaseRegionName(prop.Region))
prop=sw.prop_pTX(100E5, 100+273.15, 0.5)
print(prop.H,sw.getPhaseRegionName(prop.Region))