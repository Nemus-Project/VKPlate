%%---- test VK coefficient calc
% What is This?
% and what is it meant to do?
% who wrote it?
clear all
%close all
clc

addpath ./magpie
addpath ./vk_functions
%% ------------------------------------------------------------------------
% custom params
rho     = 8000 ;
E       = 2e11 ;
nu      = 0.3 ;
Lz      = 1e-3 ;
Lx      = 0.6 ;
Ly      = 3*0.6 ;
Nmodes  = 20 ;
npts=10;
%hvec=[0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003,0.002];
Nvec=floor(logspace(2.0,2.8,npts));
tmagpie2=zeros(npts,1);
tvk=zeros(npts,1);
%Nvec=floor(logspace(2,3.2,10));
%Nvec=[1000:100:2000];
%Hconv=zeros(Nmodes,Nmodes,Nmodes,length(Nvec));


%hfrac   = 0.005;   %-- computed as a fraction of sqrt(Lx*Ly)

%BCs Transv
BCsPhi  = [1e15 0e15 ; 1e15 0e15 ; 1e15 0e15 ; 1e15 0e15] ;
%BCsPhi  = [0 0 ; 0 0 ; 0 0 ; 0 0] ;

%BCs Airy
BCsPsi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%-- NB: these represent mathematically "clamped" BCs, but "free" physically
%-- this choice enables the "triple self-adjointness" property, so best to keep this as is

Ntensor = Nmodes;
ldim    = [Lx Ly Lz] ;
%% Reference iteration


iter=1;
h=Lx/Nvec(iter);
hvec(iter)=h;


Nxy1=Nvec(iter);

[Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,ldim,h,BCsPhi,Nmodes,"none",true) ;

Phiref=Phi;
Nxref=Nx;
Nyref=Ny;



Nxy2=Nx*Ny;

[Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,ldim,h,BCsPsi,Nmodes,"none",true) ;
Psiref=Psi;



%assert((Nxy1-Nxy2 ==0),'Nxy pure man facked')

zeta = (zetafourth).^(1/4) ;

Hv = zeros(Ntensor,Ntensor,Ntensor) ;
Ev = zeros(Ntensor,Ntensor,Ntensor) ;

Dxx = DxxBuild(Nx,Ny,h) ;
Dyy = DyyBuild(Nx,Ny,h) ;
Dxy = DxyBuild(Nx,Ny,h) ;
disp(iter)
tic;
for k = 1 : Ntensor
    Phik = Phi(:,k) ; Psik = Psi(:,k) ;

    for p = 1 : Ntensor
        Phip = Phi(:,p);

        for q = 1 : Ntensor
            Phiq = Phi(:,q) ; Psiq = Psi(:,q);
            LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;



            Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);
            %Ev(k,p,q) = trapzIntcalc(Phik.*LPhipPsiq,h,Nx,Ny)/sqrt(Phiknorm*Phipnorm*Psiqnorm) ;




        end
    end
end
toc


%%
ntest=1;
ktest=ntest;
ptest=4;
qtest=8;
ite1=2;
ite2=1;
figure
for iter=2:npts
    %------------------------------------------------------------------------

    %tic
    %------------------------------------------------------------------------
    % derived params
    %h       = sqrt(Lx*Ly)*hfrac ;
    %h       = sqrt(Lx*Ly)*hvec(iter) ;
    h=Lx/Nvec(iter);
    hvec(iter)=h;


    Nxy1=Nvec(iter);

    [Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,ldim,h,BCsPhi,Nmodes,"none",true) ;
        Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly);
    
    Nxy2=Nx*Ny;

    [Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,ldim,h,BCsPsi,Nmodes,"none",true) ;
   
     Psi = eigensign(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly);
    

    disp(iter)
    %assert((Nxy1-Nxy2 ==0),'Nxy pure man facked')

    zeta = (zetafourth).^(1/4) ;

    Hv = zeros(Ntensor,Ntensor,Ntensor) ;
    Ev = zeros(Ntensor,Ntensor,Ntensor) ;

    Dxx = DxxBuild(Nx,Ny,h) ;
    Dyy = DyyBuild(Nx,Ny,h) ;
    Dxy = DxyBuild(Nx,Ny,h) ;

    tic;
    for k = 1 : Ntensor
        Phik = Phi(:,k) ; Psik = Psi(:,k) ;
        %Phiknorm   = trapzIntcalc(Phik.*Phik,h,Nx,Ny);
        %Psiknorm   = trapzIntcalc(Psik.*Psik,h,Nx,Ny);
        for p = 1 : Ntensor
            Phip = Phi(:,p);
            %Phipnorm   = trapzIntcalc(Phip.*Phip,h,Nx,Ny) ;
            for q = 1 : Ntensor
                Phiq = Phi(:,q) ; Psiq = Psi(:,q);

                %Phiqnorm   = trapzIntcalc(Phiq.*Phiq,h,Nx,Ny) ;
                %Psiqnorm   = trapzIntcalc(Psiq.*Psiq,h,Nx,Ny) ;

                LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;

                %LPhipPsiq = vkOperator(Phip,Psiq,Dxy,Dxx,Dyy) ;

                Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);%/sqrt(Psiknorm*Phipnorm*Phiqnorm) ;
                %Ev(k,p,q) = trapzIntcalc(Phik.*LPhipPsiq,h,Nx,Ny)/sqrt(Phiknorm*Phipnorm*Psiqnorm) ;

                if iter == ite1

                    if k==ktest
                        if p==ptest
                            if q==qtest

                                Hmat=Psik.*LPhipPhiq;
                                %Hmat=Phiq;
                                mdShape = reshape(Hmat,[(Ny+1),(Nx+1)]) ;
                                xax = (0:Nx)*h ;
                                yax = (0:Ny)*h ;
                                [X,Y] = meshgrid(xax,yax) ;


                                subplot(3,3,iter-1)
                                pcolor(X,Y,3000*(mdShape));%,(abs(mdShape)));%,'FaceColor','texturemap') ;
                                shading interp
                              
                            end
                        end
                    end


                end
            end
        end
    end
    tvk=toc;
    % mdShape1 = reshape(Hmat(2,2,1,:),[(Ny+1),(Nx+1)]) ;
    % clear Hmat
    %              xax = (0:Nx)*h ;
    %             yax = (0:Ny)*h ;
    %             [X,Y] = meshgrid(xax,yax) ;
    %
    %             figure
    %
    %             mesh(X,Y,3000*(mdShape1),(abs(mdShape1)),'FaceColor','texturemap') ;
    %Emat(k,p,q,:)=Phik.*LPhipPsiq;

    H111(iter)=Hv(ntest,1,1);
    H121(iter)=Hv(ntest,2,1);
    H131(iter)=Hv(ntest,3,1);
    H141(iter)=Hv(ntest,4,1);
    H211(iter)=Hv(ntest,1,2);
    H221(iter)=Hv(ntest,2,2);
    H231(iter)=Hv(ntest,3,2);
    H241(iter)=Hv(ntest,4,2);
    H311(iter)=Hv(ntest,1,3);
    H321(iter)=Hv(ntest,2,3);
    H331(iter)=Hv(ntest,3,3);
    H341(iter)=Hv(ntest,4,3);
    H411(iter)=Hv(ntest,1,4);
    H421(iter)=Hv(ntest,2,4);
    H431(iter)=Hv(ntest,3,4);
    H441(iter)=Hv(ntest,4,4);

    sumH(1:Nmodes,iter)=sum(sum(Hv(1:Nmodes,:,:)));

    Hsp=Hv;
    Hsp(abs(Hv)<1e-1)=0;
    switch iter
        case 10
            Hsp10=Hsp;
        case 9
            Hsp9=Hsp;
        case 8
            Hsp8=Hsp;
        case 7
            Hsp7=Hsp;

        case 6
            Hsp6=Hsp;


    end

    %clear Hv
    %toc
    %nvec(iter)=Nx*Ny;

end
%%
figure
subplot(4,4,1)
semilogx(hvec,(H111),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{1,1}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,2)
semilogx(hvec,(H211),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{2,1}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,3)
semilogx(hvec,(H311),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{3,1}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,4)
semilogx(hvec,(H411),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{4,1}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)


subplot(4,4,5)
semilogx(hvec,(H121),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{1,2}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,6)
semilogx(hvec,(H221),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{2,2}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,7)
semilogx(hvec,(H321),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{3,2}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,8)
semilogx(hvec,(H421),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{4,2}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,9)
semilogx(hvec,(H131),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{1,3}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,10)
semilogx(hvec,(H231),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{2,3}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,11)
semilogx(hvec,(H331),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{3,3}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)


subplot(4,4,12)
semilogx(hvec,(H431),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{4,3}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)


subplot(4,4,13)
semilogx(hvec,(H141),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{1,4}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,14)
semilogx(hvec,(H241),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{2,4}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(4,4,15)
semilogx(hvec,(H341),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{3,4}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)


subplot(4,4,16)
semilogx(hvec,(H441),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of H^' num2str(ntest) '_{4,4}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)



%%
figure
subplot(2,2,1)
semilogx(hvec,sumH(1,:),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^1_{i,j}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(2,2,2)
semilogx(hvec,sumH(2,:),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^2_{i,j}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(2,2,3)
semilogx(hvec,sumH(3,:),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^3_{i,j}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

subplot(2,2,4)
semilogx(hvec,sumH(4,:),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^4_{i,j}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)

%%
Htest=sparse(squeeze(Hsp9(20,:,:)))
figure
spy(Htest)
%%

test=Hsp10-Hsp9;
test(abs(test)<1e2)=0;

test=sparse(squeeze(test(1,:,:)))
figure
spy(test)
%%