%%---- test VK coefficient calc
% What is This?
% and what is it meant to do?
% who wrote it?
clear all
%close all
clc


%% ------------------------------------------------------------------------
% custom params
rho     = 8000 ;
E       = 2e11 ;
nu      = 0.3 ;
Lz      = 1e-3 ;
Lx      = 0.6 ;
Ly      = 0.6 ;
T       = 0 ;
Nmodes  = 50 ;
npts=2;
%hvec=[0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003,0.002];
Nvec=floor(logspace(1,2,npts));

%Nvec=floor(logspace(2,3.2,10));
%Nvec=[1000:100:2000];
%Hconv=zeros(Nmodes,Nmodes,Nmodes,length(Nvec));


%hfrac   = 0.005;   %-- computed as a fraction of sqrt(Lx*Ly)

%BCs Transv
BCsPhi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%BCsPhi  = [0 0 ; 0 0 ; 0 0 ; 0 0] ;

%BCs Airy
BCsPsi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%-- NB: these represent mathematically "clamped" BCs, but "free" physically
%-- this choice enables the "triple self-adjointness" property, so best to keep this as is

Ntensor = Nmodes;
ldim    = [Lx Ly Lz] ;
%%
%parfor iter=1:npts
iter=npts;
%------------------------------------------------------------------------

%tic
%------------------------------------------------------------------------
% derived params
%h       = sqrt(Lx*Ly)*hfrac ;
%h       = sqrt(Lx*Ly)*hvec(iter) ;
h=Lx/Nvec(iter);
hvec(iter)=h;


Nxy1=Nvec(iter)^2;

[Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,T,ldim,h,BCsPhi,Nmodes,"none",true) ;


Nxy2=Nx*Ny;
tic;
[Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,T,ldim,h,BCsPsi,Nmodes,"none",true) ;

disp(iter)
assert((Nxy1-Nxy2 ==0),'Nxy pure man facked')

zeta = (zetafourth).^(1/4) ;

Hv = zeros(Ntensor,Ntensor,Ntensor) ;
Ev = zeros(Ntensor,Ntensor,Ntensor) ;

Dxx = DxxBuild(Nx,Ny,h) ;
Dyy = DyyBuild(Nx,Ny,h) ;
Dxy = DxyBuild(Nx,Ny,h) ;


for k = 1 : Ntensor
    Phik = Phi(:,k) ; Psik = Psi(:,k) ;
    %Phiknorm   = trapzIntcalc(Phik.*Phik,h,Nx,Ny);
    Psiknorm   = trapzIntcalc(Psik.*Psik,h,Nx,Ny);
    for p = 1 : Ntensor
        Phip = Phi(:,p);
        Phipnorm   = trapzIntcalc(Phip.*Phip,h,Nx,Ny) ;
        for q = 1 : Ntensor
            Phiq = Phi(:,q) ; Psiq = Psi(:,q);

            Phiqnorm   = trapzIntcalc(Phiq.*Phiq,h,Nx,Ny) ;
            %Psiqnorm   = trapzIntcalc(Psiq.*Psiq,h,Nx,Ny) ;

            LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;
            %LPhipPsiq = vkOperator(Phip,Psiq,Dxy,Dxx,Dyy) ;

            Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny)/sqrt(Psiknorm*Phipnorm*Phiqnorm) ;
            %Ev(k,p,q) = trapzIntcalc(Phik.*LPhipPsiq,h,Nx,Ny)/sqrt(Phiknorm*Phipnorm*Psiqnorm) ;
            
            %Hmat(k,p,q,:)=Psik.*LPhipPhiq;
            
        end
    end
end

%end
%%
% path = '/Users/alexis/Documents/MATLAB/VKPlate/src/param/';
% filename = 'NLParameters10modes.mat';
% 
% save([path filename],'rho','E','nu','Lx','Ly','Lz','Nmodes','Nx','Ny','h','BCsPhi','BCsPsi','Om','Phi','Om2','Psi','zeta','Hv')
%%







