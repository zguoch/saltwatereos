load H2O.so H2O
puts "Test tcl API"
cH2O water
set p_boling [water P_Boiling 200]
puts $p_boling