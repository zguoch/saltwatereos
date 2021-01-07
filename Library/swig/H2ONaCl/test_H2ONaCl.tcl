load H2ONaCl H2ONaCl
puts "Test tcl API"
cH2ONaCl sw
set rho [sw Rho_brine 200 300 0.2]
puts $rho
set prop [sw prop_pTX 200E5 473.15 0.2]
puts [$prop cget -H]