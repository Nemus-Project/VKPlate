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
E       = 2e11 ;
rho     = 8000 ; 
nu      = 0.3 ;
Lz      = 1e-3 ;
Lx      = 0.6 ;
Ly      = 2*0.6 ;

% Numerical parameters
Nmodes  = 10 ;
Nx=500;

% Derived values
Nvec=[100 Nx];
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

    [Om,Phi,Nx,Ny,~,~] = magpie(rho,E,nu,ldim,h,BCsPhi,Nmodes,"none",true) ; % see magpie doc

    if iter ==1

        Phiref=Phi; % Setting a refference for Psi
        Nxref=Nx; % Setting a refference for Nx
        Nyref=Ny; % Setting a refference for Ny

    else

        [Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om); % Ensures that mode order is consistent

        Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly); % Ensures that the polarization is consistent

    end

    [Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,ldim,h,BCsPsi,Nmodes,"none",false) ;% see magpie doc

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

    Dxx = DxxBuild(Nx,Ny,h) ;
    Dyy = DyyBuild(Nx,Ny,h) ;
    Dxy = DxyBuild(Nx,Ny,h) ;

    tic;
    for k = 1 : Nmodes

        Phik = Phi(:,k) ; Psik = Psi(:,k) ;

        for p = 1 : Nmodes
            
            Phip = Phi(:,p);

            for q = 1 : Nmodes

                Phiq = Phi(:,q) ; Psiq = Psi(:,q);

                LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;

                Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny); %Coupling coefficient tensor
            end
        end
    end
end
%% Defining the numerical space
xax = (0:Nx)*h ;
yax = (0:Ny)*h ;
[X,Y] = meshgrid(xax,yax) ;

%% Save parameters

if ~exist("./param/", 'dir')
    mkdir("./param/")
end

filename='parameter_clamped_1';

save(['./param/' filename '.mat'],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv');
