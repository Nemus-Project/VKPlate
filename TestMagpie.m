clear all
close all
clc
%%
% This code computes the modes of the displacement and Airy's stress of a
% nonlinear Föpple-von Karmàn plate, and the coupling coefficients between
% the modes. It takes as an input:
%
%   PHYSICAL PARAMETERS
%       E - Young modulus                           | float
%       rho - Density                               | float
%       nu - Poisson's ratio                        | float
%       Lx - x length of the plate                  | float
%       Ly - y length of the plate                  | float
%       Lz - Thickness of the plate                 | float
%       BCsPhi - Displacement boundary conditions   | 4 x 2 vector
%       BCsPsi - Stress boundary conditions         | 4 x 2 vector
%
%   NUMERICAL PARAMETERS
%       
%       Nx - Number of points in the x direction    | float
%       Nvec - Vector containing a test Nx and Nx   | 1 x 2 vector
%       Nmodes  - Number of modes                   | float
%
% The code then uses the magpie submodule to compute the linear modes of
% the plate for both displacement and stress, compare them with a
% refference through the functions, eigenmac and eigensign, and compute the
% nonlinear coupling coefficient through vkoperator.
%
% It outputs :
%
%       Hv - Coupling coeficients                   | Nmodes x Nmodes x Nmodes matrix
%       Om - Angular frequency of the Phi modes     | Nmodes x 1 vector
%       Om2 - Angular frequency of the Psi modes    | Nmodes x 1 vector
%       zetafourth - Eigenvalues of the stress      | Nmodes x 1 vector
%       Phi - Eigenvectors of displacement          | (Nx+1*Ny+1) x Nmodes matrix
%       Psi - Eigenvectors of stress                | (Nx+1*Ny+1) x Nmodes matrix
%
% The last section saves all the important variables in a file located in
% the 'param' folder, for later use by the 'main.m' code.


%% Adding path to the submodule magpie

addpath ./private/magpie

%% Variable declaration

% Physical parameters


E       = 4e+9 ;
rho     = 900 ; 
nu      = 0.4 ;
Lz      = 8e-5 ;
Lx      = 4.3e-2 ;
Ly      = Lx;
D       = E * Lz^3 / 12 / (1-nu^2);
sigma   = rho*Lz;
T       = 0%2;%2;


% E       =  200000000000.0 ;
% rho     = 8000 ; 
% nu      = 0.3 ;
% Lz      = 0.001 ;
% Lx      =  0.2 ;
% Ly      = 0.4 ;
% D       = E * Lz^3 / 12 / (1-nu^2);
% T       = 0*D;

% Numerical parameters
Nmodes  = 60 ;
Nx=300;

% Derived values
Nvec=[200 Nx];
npts=length(Nvec);


% BCs Displacement
BCsPhi  = [1e15 0e15 ; 1e15 0e15 ; 1e15 0e15 ; 1e15 0e15] ;

%BCs Airy's stress
BCsPsi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%-- NB: these represent mathematically "clamped" BCs, but "free" physically
%-- this choice enables the "triple self-adjointness" property, so best to keep this as is

ldim    = [Lx Ly Lz] ;
%%
for iter=1:npts

    h=Lx/Nvec(iter); % Defining h

    [Om,Phi,Nx,Ny,~,~] = magpie(rho,E,nu,T,ldim,h,BCsPhi,Nmodes,"none",true) ; % see magpie doc

    if iter ==1

        Phiref=Phi; % Setting a refference for Psi
        Nxref=Nx; % Setting a refference for Nx
        Nyref=Ny; % Setting a refference for Ny

    else

        [Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om); % Ensures that mode order is consistent

        Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly); % Ensures that the polarization is consistent

    end

    [Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,T,ldim,h,BCsPsi,Nmodes,"none",true) ;% see magpie doc

    if iter ==1

        Psiref=Psi; % Setting a refference for Psi

    else

        [Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly,Om2);% Ensures that mode order is consistent

        Psi = eigensign(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly);% Ensures that the polarization is consistent

    end

    disp(iter)

    zeta = (zetafourth).^(1/4) ;

    Hv = zeros(Nmodes,Nmodes,Nmodes) ;

    Ev = zeros(Nmodes,Nmodes,Nmodes) ;

    tic;
   
end
%% Defining the numerical space
xax = (0:Nx)*h ;
yax = (0:Ny)*h ;
[X,Y] = meshgrid(xax,yax) ;

%%
clear Omana OMana
mmax=20;
nmax=20;
for m=1:mmax
    for n=1:nmax
    OMana(m,n)=sqrt(D/(rho*Lz))*((m*pi/Lx)^2+(n*pi/Ly)^2);
    end
end
Omana(:)=reshape(OMana,[mmax*nmax,1]);

Omana=sort(Omana);


%% Figure
Omplt=linspace(1,60,Nmodes);

figure
stem(Omplt,Omana(1:Nmodes)/2*pi,"LineWidth",4,"MarkerSize",10,"Marker","x")
hold on
plot(Omplt,Om/2*pi,"Marker","o","MarkerSize",15,"LineWidth",4)
xlabel("Mode number")
ylabel("Frequency (Hz)")
legend("Analytical","Numerical")

set(gca,"FontSize",30)


%%
ero=abs(Omana(1:Nmodes)'-Om(:))./(2*pi*Omana(1:Nmodes)');

figure
plot(Omplt,ero,"Marker","o","MarkerSize",15,"LineWidth",4)
xlabel("Mode number")
ylabel("Relative error")
set(gca,"FontSize",30)

