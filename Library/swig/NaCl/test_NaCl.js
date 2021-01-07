var NaCl = require("./build/Release/NaCl");

var cNaCl = NaCl.cNaCl;
var T_Triple = NaCl.T_Triple;
console.log(T_Triple)
var salt = new cNaCl();
var p_boiling = salt.P_Boiling(600)
console.log("Boiling pressure at 600 deg.C is ", p_boiling)
