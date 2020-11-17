addpath('water_viscosity');
clc;clear;

addpath('water_density');

P=0.1:0.1:10;
T=0:1:300;
P=P*1E5;
[PP,TT]=meshgrid(T,P); 

PROP_SP=cal_Spivey(P,T,0);
figure(102); clf;
pcolor(PP,TT, reshape(PROP_SP.RHO(:,1,:), [length(P), length(T)])); 
shading interp; title('SP: RRHO'); colorbar;

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