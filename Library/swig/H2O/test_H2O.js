var H2O = require("./build/Release/H2O");

var cH2O = H2O.cH2O;
var aa = H2O.R_const;
console.log(aa)
var water = new cH2O();
var p_boiling = water.P_Boiling(200)
console.log("Boiling pressure at 200 deg.C is ", p_boiling)
