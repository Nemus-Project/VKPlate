clear all
close all
clc
%% ------------------------------------------------------------------------
% custom params
rho     = 7860 ;
E       = 2e+11 ;
nu      = 0.3 ;
Lz      = 4e-3 ;
Lx      = 40e-2 ;
Ly      = 80e-2 ;
T       = 0 ;

Nmodes  =10;%number of modes to compute
%With the given parameters, the case npts=10 should work, but npts=15 won't
npts=40; %number of grid size tested

Nvec=floor(logspace(1.5,2.2,npts));
tmagpie2=zeros(npts,1);
tvk=zeros(npts,1);


%BCs Transv
BCsPhi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;


%BCs Airy
BCsPsi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%-- NB: these represent mathematically "clamped" BCs, but "free" physically
%-- this choice enables the "triple self-adjointness" property, so best to keep this as is

Ntensor = Nmodes;
ldim    = [Lx Ly Lz] ;
%%
nxtestvec=ones(npts,1);
ntest=2;
ktest=ntest;
ptest=2;
qtest=3;
ltest=1;
mtest=1;
ite1=2;
ite2=1;

for iter=1:npts
    %------------------------------------------------------------------------

    %tic
    %------------------------------------------------------------------------
    % derived params
    h=Lx/Nvec(iter);
    hvec(iter)=h;
    
    


    Nxy1=Nvec(iter);

    [Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,T,ldim,h,BCsPhi,Nmodes,"none",true) ;%calls magpie

  
    
    [Phiort,OrtN]=GramSchmidt(Phi); %orthogonalise the basis of phi

    for nQ = 1 : Nmodes % normalize the basis
        Phitemp   = Phiort(:,nQ) ;
        Phinorm   = trapzIntcalc(Phitemp.*Phitemp,h,Nx,Ny) ;
        Phitemp   = Phitemp / sqrt(Phinorm) ;
        Phi(:,nQ) = Phitemp ;
    end
    clear Phitemp



    if iter ==1
        Phiref=Phi;
        Nxref=Nx;
        Nyref=Ny;

    else

        [Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om);

       

        for nQ = 1 : Nmodes % normalize the basis
            Phitemp   = Phi(:,nQ) ;
            Phinorm   = trapzIntcalc(Phitemp.*Phitemp,h,Nx,Ny) ;
            Phitemp   = Phitemp / sqrt(Phinorm) ;
            Phi(:,nQ) = Phitemp ;
        end
        clear Phitemp



        Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,T,Nmodes,Lx,Ly);


    end

    Nxy2=Nx*Ny;

    [Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,T,ldim,h,BCsPsi,Nmodes,"none",true) ;

    [Psiort,OrtN]=GramSchmidt(Psi); %orthogonalise the basis of phi

    for nQ = 1 : Nmodes % normalize the basis
        Psitemp   = Psiort(:,nQ) ;
        Psinorm   = trapzIntcalc(Psitemp.*Psitemp,h,Nx,Ny) ;
        Psitemp   = Psitemp / sqrt(Psinorm) ;
        Psi(:,nQ) = Psitemp ;
    end
    clear Phitemp


    if iter ==1
        Psiref=Psi;
    else

        [Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly,Om2);


        
        
        for nQ = 1 : Nmodes % normalize the basis
            Psitemp   = Psi(:,nQ) ;
            Psinorm   = trapzIntcalc(Psitemp.*Psitemp,h,Nx,Ny) ;
            Psitemp   = Psitemp / sqrt(Psinorm) ;
            Psi(:,nQ) = Psitemp ;
        end
        clear Phitemp

       

        Psi = eigensign(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly);

       
    end

    disp(iter)
    

    zeta = (zetafourth).^(1/4) ;

    Hv = zeros(Ntensor,Ntensor,Ntensor) ;
    

    Dxx = DxxBuild(Nx,Ny,h) ;
    Dyy = DyyBuild(Nx,Ny,h) ;
    Dxy = DxyBuild(Nx,Ny,h) ;

    tic;
    for k = 1 : Ntensor
        Phik = Phi(:,k) ; Psik = Psi(:,k) ;
       
        for p = 1 : Ntensor
            Phip = Phi(:,2);
            Phipnorm   = trapzIntcalc(Phip.*Phip,h,Nx,Ny) ;
            for q = p : Ntensor
                Phiq = Phi(:,q) ; Psiq = Psi(:,q);

               
                LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;

                
                Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);
                
                Hv(k,q,p) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);
                
            end
        end
    end
   


   
  if Nx*h-Nxref*hvec(1) ~= 0
        disp("This point is fudged")
        nxtestvec(iter)=0;
      
    end
    clear Hv
   

   

end


