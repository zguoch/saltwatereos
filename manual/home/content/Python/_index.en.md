---
title: "Python"
date: 2018-12-29T11:02:05+06:00
icon: "ti-panel"
logo: "images/python.svg"
# description: ""
type : "api"
---

## Install [pyswEOS](https://pypi.org/project/pyswEOS/)

<!-- <div>
<select id="pyswEOS-os" class="form-select" aria-label="pyswEOS-OS" onchange="pyswEOSselect(this);">
  <option selected>Select OS</option>
  <option value="1">MacOSX</option>
  <option value="2">Linux</option>
  <option value="3">Windows</option>
</select>
<select id="pyswEOS-pyversion" class="form-select" aria-label="pyswEOS-pyversion" onchange="pyswEOSselect(this);">
  <option selected>Select Python Version</option>
  <option value="1" id="dadfs">2.7</option>
  <option value="2">3.5</option>
  <option value="3">3.6</option>
  <option value="3">3.7</option>
  <option value="3">3.8</option>
  <option value="3">3.9</option>
</select>
<div class="btn-group">
  <button class="btn btn-secondary btn-sm" type="button">
    Small split button
  </button>
  <button type="button" class="btn btn-sm btn-secondary dropdown-toggle dropdown-toggle-split" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
    <span class="sr-only">Toggle Dropdown</span>
  </button>
  <div class="dropdown-menu" id="pyswEOS-dropdown">
    <a class="dropdown-item" href="#">Action</a>
    <a class="dropdown-item" href="#">Action</a>
  </div>
</div>
</div> -->

```bash
pip install pyswEOS
```

## How to use: demo

<h4 class="text-center">Please visit <a href="../manual/index.html">online manual</a> for more details.</h4>

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