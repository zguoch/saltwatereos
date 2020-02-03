%% function that computes density of pure water as in Spivey et al., 2004
% equations (4)-(9)


function [ rho_w ] = density_pure_water(t,p)

p0      = 70e6;                                                     % [Pa] Reference pressure

%% Get coefficients for pure water
[ rho_w0, Ew, Fw ]   = coefficients_pure_water( t );

%% Compute isothermal compressibility
 Iw                  = isothermal_compressibility(p,p0,Ew,Fw);
 Iw0                 = isothermal_compressibility(p0,p0,Ew,Fw);

%% Compute density of pure water
 rho_w               = rho_w0 .* exp(Iw - Iw0);                       % eq. (8) Spivey et al., 2004
rho_w                = rho_w.*1000;
end

%%--------------------------------------------------------------------------
%%                             References
%%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation volume factor, compressibility,
% methane solubility, and viscosity for oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology, 43 (7), 52-60.