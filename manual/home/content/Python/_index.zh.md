---
title: "Python"
date: 2018-12-29T11:02:05+06:00
icon: "ti-pie-chart"
logo: "images/python.png"
# description: ""
type : "api"
---


## 安装[pyswEOS](https://pypi.org/project/pyswEOS/)

```bash
pip install pyswEOS
```

## 使用示例

<h4 class="text-center">请访问<a href="../manual/index.html">在线文档</a>获取更多细节信息</h4>

* H2O

```py
from pyswEOS import H2O
water=H2O.cH2O()
print(water.P_Boiling(200))
```

* NaCl

```py
from pyswEOS import NaCl
salt=NaCl.cNaCl()
print('Triple temperature: ',NaCl.T_Triple)
print(salt.P_Boiling(900))
```

* H2ONaCl

```py
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
```