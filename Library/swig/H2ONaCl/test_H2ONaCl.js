var H2ONaCl = require("./build/Release/H2ONaCl");

var cH2ONaCl = H2ONaCl.cH2ONaCl;
var TMAX_K = H2ONaCl.TMAX_K;
console.log(TMAX_K)
var sw = new cH2ONaCl();
var x_wt = sw.Mol2Wt(0.1);
console.log("0.1 mol/kg is  ", x_wt, " wt% NaCl")
var x_mol = sw.Wt2Mol(x_wt)
console.log(x_wt, "wt% NaCl is  ", x_mol, " mol/kg")
