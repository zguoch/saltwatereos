load H2ONaCl H2ONaCl
puts "Phase Region"
cH2ONaCl sw
set region [sw findPhaseRegion 400 200 0.032]
puts $region
set region_name [sw getPhaseRegionName $region]
puts $region_name