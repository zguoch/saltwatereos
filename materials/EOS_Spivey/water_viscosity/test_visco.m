% Test viscosity brine 

T_ = 20+273.15:1:275+273.15;
p  = 100e6; 
c  = 0; 


mu_ = zeros(numel(T_),1); 


for i = 1:numel(T_)
  t = T_(i);  
  [ mu_b ] = viscosity_brine( p,t,c );
  mu_(i) = mu_b;
end 
figure()
plot(mu_.*1e6)