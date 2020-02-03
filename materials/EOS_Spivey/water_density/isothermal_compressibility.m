%% Fonction that computes the isothermal compressibility
%% equations (7) and (14) in Spivey et al., 2004

function [ I ] = isothermal_compressibility(p,p0,E,F)

I       = 1./E .* log(E.*(p./p0) + F);

end

