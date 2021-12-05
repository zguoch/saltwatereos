import numpy as np 
import random
import time
from pyswEOS import H2ONaCl
from pyswEOS import H2O
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

# 6. test lookup table
# 6.1 create lookup table
dim =2
xmin, ymin = 1 + 273.15, 5E5 #T [K], P[Pa]
xmax, ymax = 700 + 273.15, 400E5
X_wt = 0.2 #wt% NaCl [0,1]
min_level = 4
max_level = 10
# sw.createLUT_2D_PTX("constX", xmin,xmax, ymin,ymax, X_wt, min_level, max_level)
# sw.save_to_vtk("test.vtu")
# sw.save_to_binary("test.bin")

# 6.2 load from bin file and save to vtu
sw.loadLUT_PTX("lut_TPX_13.bin")
# sw.save_to_vtk("test_out.vtu")

# 6.3 lookup 
start = time.process_time()
ind = 0
for i in range(0, int(1E6)):
    T_K     = (random.random())*(xmax - xmin) + xmin
    p_Pa    = (random.random())*(ymax - ymin) + ymin
    # print(T_K, p_Pa)
    prop    = sw.searchLUT_2D_PTX(T_K, p_Pa)
    prop2   = sw.prop_pTX(p_Pa, T_K, X_wt)
    # print(prop.Region, prop2.Region)
    if(prop.Region != prop2.Region):
        # print('diff %d: %d'%(ind, prop.Region))
        ind += 1
    else:
        err = prop.Rho - prop2.Rho
        if(np.abs(err)>0.5): print('Rho err: %f'%(err))
stop = time.process_time()
print("Search, elapsed time: %.2f, %d need refine" % (stop - start, ind))