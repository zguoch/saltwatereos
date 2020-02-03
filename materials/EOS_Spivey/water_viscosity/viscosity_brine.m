function [ mu_b ] = viscosity_brine( p,t,c )
%% function that computes density of a brine as in Spivey et al., 2004
% equations (35)
% Author: M.Collignon, Uni. Geneva, 2018
% Input parameters  
% p - pressure [Pa] 
% t - temperature [degree] 
% c - concentration [kg/kg]

% Output parameters
% rho_b - brine density  [kg/m^3]
p               = p./1e6;                                                   % conversion from Pa tp MPa 
%t               = t-273.15;                                                 % conversion from Kelvin to �Celsius
MW_NaCl         = 0.0584428;                                                % [kg/mol] Molar weigth of NaCl
X_NaCl          = c;                                                        % mass fraction of NaCl  
X_Water         = 1 - X_NaCl;                                               % mass fraction of water    
cNaCl           = X_NaCl./X_Water;                                          % concentration of NaCl per kg of water [kg/kg]   
mol_NaCl        = cNaCl./MW_NaCl;                                           % molality [mol/kg H2O]  

tMin            = 0;
tMax            = 275;
pMin            = 0;
pMax            = 200;
mol_NaClMin     = 0;
mol_NaClMax     = 5.7;

if max(double(t)) > 275
    warning('T out of range')
end

if max(double(p)) > 200
    warning('p out of range')
end

if max(double(mol_NaClMax)) > 5.7;
    warning('mol_NaCl out of range')
end

%% Caping for the ADI variable scheme
t = min(max(t, tMin), tMax);
p = min(max(p, pMin), pMax);
mol_NaCl = min(max(mol_NaCl, mol_NaClMin), mol_NaClMax);


%% Coefficient for the scaling
A1 = 0.0173; 
A2 = 0.068;
B1 = -1.0531; 
B2 = 0.0273;


mu_kestin = kestin_brine_viscosity(p, 125,mol_NaCl); % function kestin_brine_viscosity take MPa in input
mu_k      = kestin_brine_viscosity(p,t,mol_NaCl);

term1 = (A2.*(p./100)+A1).*(log(t./125)).^2;
term2 = (B2.*(p./100)+B1).*(log(t./125));

mu    = mu_kestin.*exp(term1 + term2);

if double(t) > 125 
    mu_b = mu;
else 
    mu_b = mu_k;
end
    mu_b = mu_b.*1e-6; 
end

