% compute the brine density (Spivey et al., 2004)

T = 10;                     % temperature in degrees
p = 100;                    % pressure in MPa
c = 0;                      % NaCl concentration in 

T_ = 20:1:350;

%% Coefficient for brine 

mu_ = zeros(numel(T_),1); 
mu_k_ = zeros(numel(T_),1);

for i = 1:numel(T_)
    T = T_(i);
A1 = 0.0173; 
A2 = 0.068;
B1 = -1.0531; 
B2 = 0.0273; 

mu_kestin = kestin_brine_viscosity(p, 125,c); % function kestin_brine_viscosity take MPa in input
mu_k = kestin_brine_viscosity(p,T,c);

term1 = (A2.*(p./100)+A1).*(log(T./125)).^2;
term2 = (B2.*(p./100)+B1).*(log(T./125));

mu = mu_kestin.*exp(term1 + term2);

mu_k_(i) = mu_k;

if T>125
mu_(i)   = mu;
else
mu_(i) = mu_k_(i); 
end
end 