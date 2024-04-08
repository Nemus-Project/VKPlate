clear all
close all
clc
addpath magpie
%%
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
bc=[1e4 1e6 1e8 1e10 1e15];
n=0
figure
for i=1:5
BCs(1,:)=bc(i);



%--- derived parameters (don't change here)
rho     = 8000 ;
E       = 2e11 ;
nu      = 0.3 ;
Lz      = 1e-3 ;
Lx      = 0.6 ;
Ly      = 0.4 ;
hfrac   = 0.002;   %-- computed as a fraction of sqrt(Lx*Ly)
Nmodes  = 6  ;
h       = sqrt(Lx*Ly)*hfrac ;
%%

[Om,Q,Nx,Ny,biharm] = magpie(rho,E,nu,[Lx Ly Lz],h,BCs,Nmodes,'none',true);

colormap('parula') ;
        xax = (0:Nx)*h ;
        yax = (0:Ny)*h ;
        [X,Y] = meshgrid(xax,yax) ;

        
        for m = 2 : Nmodes
            mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
            subplot(Nmodes-1,5,m+n-1)
            mesh(X,Y,3000*(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;
        end
        n=n+5
end

