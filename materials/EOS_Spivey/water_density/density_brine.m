%% function that computes density of a brine as in Spivey et al., 2004
%% equations (15)

function [ rho_b ] = density_brine(s,t,p)

%MW_NaCl         = 0.0584428;                    % [kg/mol]      Molar weigth of NaCl
%mol_NaCl        = s*1e-3/MW_NaCl;               % [mol/kg]
p0              = 70e6;                         % [Pa]          Reference pressure
MW_NaCl         = 0.0584428;                                                % [kg/mol] Molar weigth of NaCl
X_NaCl          = s;                                                        % mass fraction of NaCl
X_Water         = 1 - X_NaCl;                                               % mass fraction of water
cNaCl           = X_NaCl./X_Water;                                          % concentration of NaCl per kg of water [kg/kg]
mol_NaCl        = cNaCl./MW_NaCl;



%% Get coefficients for pure water to compute density brine (rho_b0) and coefficients for the brine (Eb, Fb)
[ rho_w0, Ew, Fw ]   = coefficients_pure_water( t );

%% Compute the density of the brine at the reference pressure
 rho_b0              = density_brine_refP(rho_w0,t,mol_NaCl);

%% Compute the coefficients Eb and Fb for the brine
[ Eb,Fb ]            = coefficients_brine(Ew,Fw,mol_NaCl,t);

%% Compute isothermal compressibility
 Ib                  = isothermal_compressibility(p,p0,Eb,Fb);
 Ib0                 = isothermal_compressibility(p0,p0,Eb,Fb);

%% Compute density
 rho_b               = rho_b0 .* exp(Ib - Ib0);                   % Eq.(15) in Spivey et al., 2004
 rho_b               = rho_b.*1000;

end

%%--------------------------------------------------------------------------
%%                             References
%%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation volume factor, compressibility,
% methane solubility, and viscosity for oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology, 43 (7), 52-60.