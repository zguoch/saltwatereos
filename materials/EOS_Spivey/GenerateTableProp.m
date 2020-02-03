%% Script that load data from table (from Driesner) 

clc;clear;
addpath('water_viscosity');
addpath('water_density');
addpath('/Users/zguo/MyData/Research/3-CodeProject/H2O_NaCl_EOS_Driesner_2007');

% 1. load data of Drisner
% PROP=GenTable;

load('PROP_Dr.mat');

figure(101); clf;
[TT,PP]=meshgrid(PROP.T, PROP.P); 
pcolor(PP,TT, reshape(PROP.RHO(:,1,:), [length(PROP.P), length(PROP.T)])); 
shading interp; title('Driesner: RRHO'); colorbar;


% 2. load table of Falko
DATA=load('PROP_FV');
PROP_FV=DATA.PROP2;

% 3. Spivey
PROP_SP=cal_Spivey(PROP.P*1e5,PROP.T,0);
figure(102); clf;
pcolor(PP,TT, reshape(PROP_SP.RHO(:,1,:), [length(PROP.P), length(PROP.T)])); 
shading interp; title('SP: RRHO'); colorbar;
figure(200); pcolor(PP,TT, reshape(abs(PROP_SP.RHO(:,1,:) - PROP.RHO(:,1,:))./PROP.RHO(:,1,:)*100, [length(PROP.P), length(PROP.T)])); shading interp; colorbar('fontsize',20);
% calculate FV
PROP.P(1)=5;
[RHO, MU]=cal_FV(PROP.T,PROP.P*1e5,0.000001);

return;

% PXT
index_p=10;
% compar_p('RHO',index_p,PROP_SP,PROP_FV, 'SP - FV');

index_x=1;
compar_x('RHO',index_x,PROP_SP,PROP, 'SP - Dr');

function [RHO, MU]=cal_FV(T,P,X)
    for i=1:length(T)
        for j=1:length(P)
            for k=1:length(X)
                [ PROP ] = fluidprop_NaCl_PTX( P(j),T(i), X(k), 1 );
                RHO(i,j)=PROP.Rho;
                MU(i,j)=PROP.mu_l;
            end
        end 
    end
    
    [PP,TT]=meshgrid(P,T);
    figure(100);clf;
    pcolor(PP/1E5,TT,RHO); shading interp;
end


function compar_x(name_prop,index_x, PROP, PROP_FV, prop1_prop2)
prop11=getfield(PROP,name_prop);
prop22=getfield(PROP_FV, name_prop);
prop1=prop11(:,index_x,:);
prop2=prop22(:,index_x,:);
P=PROP.P;
X=PROP.X(index_x);
T=PROP.T;

res=(prop1-prop2)./prop2.*100;
[TT,PP]=meshgrid(T,P);
res2=reshape(res,[size(res,1),size(res,3)]);
pcolor(PP,TT,res2); shading interp
colorbar('fontsize',20);
title([prop1_prop2, ': ' name_prop ' at fixed salinity x=' num2str(X)], 'fontsize', 20)
xlabel('Pressure (bar)', 'fontsize', 20);
ylabel('Temperature (C)', 'fontsize', 20);
end


function compar_p(name_prop,index_p, PROP, PROP_FV, prop1_prop2)
prop11=getfield(PROP,name_prop);
prop22=getfield(PROP_FV, name_prop);
prop1=prop11(index_p,:,:);
prop2=prop22(index_p,:,:);
p=PROP.P(index_p);
X=PROP.X;
T=PROP.T;

res=abs(prop1-prop2)./prop2.*100;
[TT,XX]=meshgrid(T,X);
res2=reshape(res,[size(res,2),size(res,3)]);
pcolor(XX,TT,res2); shading interp
colorbar('fontsize',20);
title([prop1_prop2, ': ' name_prop ' at fixed pressure p=' num2str(p) ' bar'], 'fontsize', 20)
xlabel('Salinity', 'fontsize', 20);
ylabel('Temperature (C)', 'fontsize', 20);
end



function PROP_SP=cal_Spivey(P,T,X)
for i=1:length(P)
    for j=1:length(X)
        for k=1:length(T)
            s=X(j);
            t=T(k);
            p=P(i);
            PROP_SP.RHO(i,j,k)  = density_brine(s,t,p);
            PROP_SP.MU(i,j,k)=viscosity_brine(p,t,s);
        end
    end
    PROP_SP.P=P;
    PROP_SP.T=T;
    PROP_SP.X=X;
end
end

function PROP=GenTable()

%% Load data 
[data,txt] = xlsread('DriesnerTable','param');


%% Reshape data
% get the size and range for p, T, and salinity

Trange      =   5:5:150; 
Prange      =   0:25:500;
Prange(1)   =   1; 
srange      =   0:0.01:0.25;

nT          =   length(Trange);
nP          =   length(Prange);
ns          =   length(srange);
ntot        =   nT*ns*nP;

% extract param from data 
rho         =   data(:,3);
Cp          =   data(:,4);
hfl         =   data(:,5);
mu          =   data(:,6);

% reshape value
rho_mat     = reshape(rho, nP, ns, nT);
Cp_mat      = reshape(Cp, nP, ns, nT);
hfl_mat     = reshape(hfl, nP, ns, nT);
mu_mat      = reshape(mu, nP, ns, nT);

PROP.RHO=rho_mat;
PROP.CP=Cp_mat;
PROP.H=hfl_mat;
PROP.MU=mu_mat;
PROP.T=Trange;
PROP.P=Prange;
PROP.X=srange;
end
% 
% %% Get properties at fixed reference parameters and generate plot 
% % get densities
% rho_pt      = rho_mat(:,1,:); 
% rho_pt      = reshape(rho_pt, nP, nT);
% rho_ps      = rho_mat(:,:,10);
% rho_ts      = rho_mat(3,:,:); 
% rho_ts      = reshape(rho_ts, ns, nT);
% 
% % get viscosity
% mu_pt      = mu_mat(:,1,:); 
% mu_pt      = reshape(mu_pt, nP, nT);
% mu_ps      = mu_mat(:,:,10);
% mu_ts      = mu_mat(3,:,:); 
% mu_ts      = reshape(mu_ts, ns, nT);
% 
% % get heat capacity 
% Cp_pt      = Cp_mat(:,1,:); 
% Cp_pt      = reshape(Cp_pt, nP, nT);
% Cp_ps      = Cp_mat(:,:,10);
% Cp_ts      = Cp_mat(3,:,:); 
% Cp_ts      = reshape(Cp_ts, ns, nT);
% 
% % get enthalpy 
% hfl_pt      = hfl_mat(:,1,:); 
% hfl_pt      = reshape(hfl_pt, nP, nT);
% hfl_ps      = hfl_mat(:,:,10);
% hfl_ts      = hfl_mat(3,:,:); 
% hfl_ts      = reshape(hfl_ts, ns, nT);
% 
% %% get grid
% 
% [T1, p1] = meshgrid(Trange, Prange);
% [s2, p2] = meshgrid(srange, Prange); 
% [T3, s3] = meshgrid(Trange, srange); 
% 
% %% plotting
% figure()
% 
% % plot density
% sx1 = subplot(4,3,1)
% pcolor(p1, T1, rho_pt)
% shading interp 
% colorbar
% xlabel('pressure [bar]')
% ylabel('temperature [ºC]')
% title('Density [kg m^{-3}] , X_{NaCl} = 0')
% 
% sx2 = subplot(4,3,2)
% pcolor(p2, s2, rho_ps)
% shading interp 
% colorbar
% xlabel('pressure [bar]')
% ylabel('X_{NaCl}')
% title('Density [kg m^{-3}] , T = 50ºC')
% 
% sx3 = subplot(4,3,3)
% pcolor(T3, s3, rho_ts)
% shading interp 
% colorbar
% xlabel('temperature [ºC]')
% ylabel('X_{NaCl}')
% title('Density [kg m^{-3}] , P = 50 bar')
% 
% 
% % plot viscosity
% sx4 = subplot(4,3,4)
% pcolor(p1, T1, mu_pt)
% shading interp 
% colorbar
% colormap(sx4,'hot')
% xlabel('pressure [bar]')
% ylabel('temperature [ºC]')
% title('Viscosity [Pa s] , X_{NaCl} = 0')
% 
% sx5 = subplot(4,3,5)
% pcolor(p2, s2, mu_ps)
% shading interp 
% colorbar
% colormap(sx5,'hot')
% xlabel('pressure [bar]')
% ylabel('X_{NaCl}')
% title('Viscosity [Pa s] , T = 50ºC')
% 
% sx6 = subplot(4,3,6)
% pcolor(T3, s3, mu_ts)
% shading interp 
% colorbar
% colormap(sx6,'hot')
% xlabel('temperature [ºC]')
% ylabel('X_{NaCl}')
% title('Viscosity [Pa s] , P = 50 bar')
% 
% % plot heat capacity 
% sx7 = subplot(4,3,7)
% pcolor(p1, T1, Cp_pt)
% shading interp 
% colorbar
% colormap(sx7,'summer')
% xlabel('pressure [bar]')
% ylabel('temperature [ºC]')
% title('Heat Capacity [J kg^{-1} K^{-1}] , X_{NaCl} = 0')
% 
% sx8 = subplot(4,3,8)
% pcolor(p2, s2, Cp_ps)
% shading interp 
% colorbar
% colormap(sx8,'summer')
% xlabel('pressure [bar]')
% ylabel('X_{NaCl}')
% title('Heat Capacity [J kg^{-1} K^{-1}] , T = 50ºC')
% 
% sx9 = subplot(4,3,9)
% pcolor(T3, s3, Cp_ts)
% shading interp 
% colorbar
% colormap(sx9,'summer')
% xlabel('temperature [ºC]')
% ylabel('X_{NaCl}')
% title('Heat Capacity [J kg^{-1} K^{-1}] , P = 50 bar')
% 
% % plot enthalpy
% sx10 = subplot(4,3,10)
% pcolor(p1, T1, hfl_pt)
% shading interp 
% colorbar
% colormap(sx10,'cool')
% xlabel('pressure [bar]')
% ylabel('temperature [ºC]')
% title('Enthalpy [J kg^{-1}] , X_{NaCl} = 0')
% 
% sx11 = subplot(4,3,11)
% pcolor(p2, s2, hfl_ps)
% shading interp 
% colorbar
% colormap(sx11,'cool')
% xlabel('pressure [bar]')
% ylabel('X_{NaCl}')
% title('Enthalpy [J kg^{-1}] , T = 50ºC')
% 
% sx12 = subplot(4,3,12)
% pcolor(T3, s3, hfl_ts)
% shading interp 
% colorbar
% colormap(sx12,'cool')
% xlabel('temperature [ºC]')
% ylabel('X_{NaCl}')
% title('Enthalpy [J kg^{-1}] , P = 50 bar')
% 
% %% load spivey
% 
% % Sp        = load('Sp_param.mat');
% Sp_mu_ps  = Sp.mu_ps;
% Sp_mu_pt  = Sp.mu_pT;
% Sp_mu_ts  = Sp.mu_ts;
% 
% Sp_rho_ps = Sp.rho_ps;
% Sp_rho_pt = Sp.rho_pT;
% Sp_rho_ts = Sp.rho_ts;
% 
% %% compute residual 
% 
% r_mu_ps   = 100.*abs(Sp_mu_ps-mu_ps)./mu_ps;
% r_mu_pt   = 100.*abs(Sp_mu_pt-mu_pt)./mu_pt;
% r_mu_ts   = 100.*abs(Sp_mu_ts-mu_ts)./mu_ts;
% 
% r_rho_ps   = 100.*abs(Sp_rho_ps-rho_ps)./rho_ps;
% r_rho_pt   = 100.*abs(Sp_rho_pt-rho_pt)./rho_pt;
% r_rho_ts   = 100.*abs(Sp_rho_ts-rho_ts)./rho_ts;
% 
% %% plot residual 
% figure()
% 
% % plot density
% sx1 = subplot(2,3,1)
% pcolor(p1, T1, r_rho_pt)
% shading interp 
% colorbar
% xlabel('pressure [bar]', 'FontSize', 14)
% ylabel('temperature [ºC]', 'FontSize', 14)
% title('Density residual [%] , X_{NaCl} = 0', 'FontSize', 14)
% set(gca,'FontSize', 12)
% 
% sx2 = subplot(2,3,2)
% pcolor(p2, s2, r_rho_ps)
% shading interp 
% colorbar
% xlabel('pressure [bar]', 'FontSize', 14)
% ylabel('X_{NaCl}', 'FontSize', 14)
% title('Density residual [%] , T = 50ºC', 'FontSize', 14)
% set(gca,'FontSize', 12)
% 
% sx3 = subplot(2,3,3)
% pcolor(T3, s3, r_rho_ts)
% shading interp 
% colorbar
% xlabel('temperature [ºC]','FontSize', 14)
% ylabel('X_{NaCl}', 'FontSize', 14)
% title('Density residual [%] , P = 50 bar', 'FontSize', 14)
% set(gca,'FontSize', 12)
% 
% % plot viscosity
% sx4 = subplot(2,3,4)
% pcolor(p1, T1, r_mu_pt)
% shading interp 
% colorbar
% colormap(sx4,'hot')
% xlabel('pressure [bar]')
% ylabel('temperature [ºC]')
% title('Viscosity residual [%] , X_{NaCl} = 0', 'FontSize', 14)
% set(gca,'FontSize', 12)
% 
% sx5 = subplot(2,3,5)
% pcolor(p2, s2, r_mu_ps)
% shading interp 
% colorbar
% colormap(sx5,'hot')
% xlabel('pressure [bar]')
% ylabel('X_{NaCl}')
% title('Viscosity residual [%] , T = 50ºC', 'FontSize', 14)
% set(gca,'FontSize', 12)
% 
% sx6 = subplot(2,3,6)
% pcolor(T3, s3, r_mu_ts)
% shading interp 
% colorbar
% colormap(sx6,'hot')
% xlabel('temperature [ºC]')
% ylabel('X_{NaCl}')
% title('Viscosity residual [%] , P = 50 bar', 'FontSize', 14)
% set(gca,'FontSize', 12)
% 
% %% Plotting: get a grid 
% [T,P] = meshgrid(Trange, Prange);
% 
% 
% %% plot properties for pure water
% rho_pw  = rho_mat(:,1,:); 
% rho_pw  = reshape(rho_pw, nP, nT);
% 
% figure()
% surf(P,T,rho_pw) 
% shading interp 
% axis tight 
% xlabel ('Pressure [bar]')
% ylabel ('Temperature [ºC]')
% zlabel('density [kg m^{-3}]')
% title('density of pure water')
% %% plot properties for water at 0.01 XNaCl 
% rho_nacl  = rho_mat(:,2,:); 
% rho_nacl  = reshape(rho_nacl, nP, nT);
% 
% figure()
% surf(P,T,rho_nacl) 
% shading interp 
% axis tight 
% xlabel ('Pressure [bar]')
% ylabel ('Temperature [ºC]')
% zlabel('density [kg m^{-3}]')
% title('density of water at 0.01 XNaCl (eq. 10g L^{-1})')
% 
% %% plot viscosity pure water 
% 
% mu_pw  = mu_mat(:,1,:); 
% mu_pw  = reshape(mu_pw, nP, nT);
% 
% figure()
% surf(P,T,mu_pw) 
% shading interp 
% axis tight 
% xlabel ('Pressure [bar]')
% ylabel ('Temperature [ºC]')
% zlabel('viscosity [Pa s]')
% title('viscosity of pure water')
% 
% %% plot viscosity for water at 0.01 XNaCl
% mu_nacl  = mu_mat(:,2,:); 
% mu_nacl  = reshape(mu_nacl, nP, nT);
% 
% figure()
% surf(P,T,mu_nacl) 
% shading interp 
% axis tight 
% xlabel ('Pressure [bar]')
% ylabel ('Temperature [ºC]')
% zlabel('viscosity [Pa s]')
% title('viscosity of water at 0.01 XNaCl (eq. 10g L^{-1})')
% 
% %% Save each parameters in a separated matrix
% 
% save('DriesnerProps','rho_mat','mu_mat','hfl_mat','Cp_mat', 'Trange', 'Prange', 'srange');




