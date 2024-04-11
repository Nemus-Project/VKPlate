clear all
close all
clc
addpath magpie
%% Biharm case stuff

%-- elastic constants around the edges
K0y     = 1e15 ;
R0y     = 0e15 ;
Kx0     = 1e15 ;
Rx0     = 0e15 ;
KLy     = 1e15 ;
RLy     = 0e15;
KxL     = 1e15 ;
RxL     = 0e15 ;

BCs=zeros(4,2);
BCs(:,1)=1e15;



%--- derived parameters (don't change here)
Lz=1e-3;
nu=0.3;
E=1e9;
D       = E * Lz^3 / 12 / (1-nu^2) ;
%k       =1/44100;
%kappa   =sqrt(D/(rho*Lz));
%h       =2*sqrt(kappa*k)
%h       = 0.67*(sqrt(Lx*Ly)*0.01);
%Nx      = floor(Lx/h) ;
%Ny      = floor(Ly/h) ;

%bhmat(BCs,[6,7],0.01,Lz,E,nu,'diag');

%%

nbpts=100;
modes=zeros(nbpts,1);


 Nx=100;
 Ny=100;
 h=0.01;
 ncalc=10;
 


 %%
nbpts=100;
eig10t=zeros(nbpts,1);
eig20t=zeros(nbpts,1);
eig50t=zeros(nbpts,1);

nxny=zeros(nbpts,1);
for n=nbpts:-1:1
    Nx=10*n;
    Ny=1000;
    h=0.01;
    nxny(n)= Nx*Ny;
    biHarm=bhmat(BCs,[Nx,Ny],h,Lz,E,nu,'diag');

    tic
    [~,~] = eigs(biHarm,10,'smallestabs') ;
    eig10t(n)=toc

    tic
    [~,~] = eigs(biHarm,20,'smallestabs') ;
    eig20t(n)=toc

    tic
    [~,~] = eigs(biHarm,50,'smallestabs') ;
    eig50t(n)=toc
    disp(n)


end
%%
figure
plot(nxny,eig10t,nxny,eig20t,nxny,eig50t,Linewidth=3)
legend("10 modes","20 modes","50 modes")
xlabel("Total Number of Grid Points")
ylabel("Time (s)")
set(gca,'Fontsize',20)