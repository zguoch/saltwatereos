import H2ONaCl
import H2O
water=H2O.cH2O()
sw=H2ONaCl.cH2ONaCl()
# 1. function
print(sw.Rho(300, 200, 0.2))
print(sw.P_VaporLiquidHaliteCoexist(200))
# 2. function with multi returns
p_crit,x_crit=sw.P_X_Critical(800)
print(p_crit,x_crit,sw.Mol2Wt(x_crit)*100)
# 3. constants
print(H2ONaCl.TMAX)
# 4. 
print(water.P_Boiling(200))