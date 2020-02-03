%% compute the viscosity of brine according to Kestin et al. 1978
% 

c = 1;                         % concentration in NaCl in mole/kg H2O. 
T_ = 20:5:150;                             % temperature in �C
p = 0.1;                           % pressure in MPa


 mu_ = zeros(numel(T_),1); 
 
for i = 1:numel(T_)
    T = T_(i);
%% ---- compute the beta factors ----
% compute reduced concentration c_star = c/cs 
cs = 6.044 + (0.28e-2).*T + (0.36e-4).*(T.^2); 
c_star = c./cs; 

% pressure coefficient for water
beta_w = -1.297 + (0.574e-1).*T - (0.697e-3).*(T.^2) + (0.447e-5).*(T.^3) - (0.105e-7).*(T.^4);

% excess pressure coefficient at the saturation 
beta_sE = 0.545 + (0.28e-2).*T - beta_w; 

% reduced coefficient for pressure
beta_star = 2.5*c_star - 2.0*(c_star.^2) + 0.5*(c_star.^3);

% beta factor for viscosity equation 
beta = beta_sE.*beta_star + beta_w;

% beta_v(i) = beta;
% end 
%% ---- compute the different viscosities ----

% Compute A and B factors 

   A = (0.3324e-1).*c + (0.3624e-2).*(c^2) - (0.1879e-3).*(c^3);
   B = -(0.396e-1).*c + (0.102e-1).*(c^2) - (0.702e-3).*(c^3);

%  A = (0.3324e-1 + 0.3624e-2 - 0.1879e-3).*(c^2);
%  B = (-0.396e-1 + 0.102e-1 - 0.702e-3).*(c^2);



mu_w0_20 = 1002;                                                           % viscosity of mu_w0 at 20�C in �Pa s. 

% compute �_w0 (mu_w0) 
log_mu_w0_ov_mu_w0_20 = (20-T)./(96+T).*[1.2378 - (1.303e-3).*(20-T) + (3.06e-6).*((20-T).^2) ...
                                        + (2.55e-8).*((20-T).^3)];

%mu_w0 = exp(log_mu_w0_ov_mu_w0_20).*mu_w0_20;   
mu_w0 = (10^log_mu_w0_ov_mu_w0_20).*mu_w0_20;

% compute �_0 (mu_0)
log_mu_0_ov_mu_W0 = A + B.*log_mu_w0_ov_mu_w0_20;

%mu_0 = exp(log_mu_0_ov_mu_W0).*mu_w0;                                    
mu_0 = (10^log_mu_0_ov_mu_W0).*mu_w0;                                  

% finally compute the viscosity of brine 
mu = mu_0.*(1+beta.*p.*1e-3);

mu_(i) = mu; 


end
