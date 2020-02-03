% test the density of the model of Spivey et al., 2004 
% compare with the data from Duan and Mao Table 15

clear all, close all 

t         = 100; 
p         = 2000e5;
MW_NaCl   = 0.0584428;
mol_NaCl  = 6; 
s         = mol_NaCl*MW_NaCl*1000; 

format long
% rho_w  = density_pure_water(t,p)
% rho_b  = density_brine(s,t,p)

Trange      =   5:5:150; 
Prange      =   0:25:500;
Prange(1)   =   1; 
srange      =   0:0.01:0.25;
[t,p,x]=meshgrid(Trange, Prange*1e5, srange);

RHO  = density_brine(s,t,p);