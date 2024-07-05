%%---- test VK coefficient calc
% What is This?
% and what is it meant to do?
% who wrote it?
clear all
%close all
clc

addpath ./private/magpie

%% ------------------------------------------------------------------------
% custom params
rho     = 8000 ;
E       = 2e11 ;
nu      = 0.3 ;
Lz      = 1e-3 ;
Lx      = 0.6 ;
Ly      = 2*0.6 ;
Nmodes  = 10 ;
npts=Nmodes;
%hvec=[0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003,0.002];
Nvec=floor(logspace(2,2.5,npts));
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
%%
ntest=1;
ktest=ntest;
ptest=5;
qtest=2;
ltest=1;
mtest=1;
ite1=2;
ite2=1;
figure
for iter=1:npts
    %------------------------------------------------------------------------

    %tic
    %------------------------------------------------------------------------
    % derived params
    %h       = sqrt(Lx*Ly)*hfrac ;
    %h       = sqrt(Lx*Ly)*hvec(iter) ;
    h=Lx/Nvec(iter);
    hvec(iter)=h;


    Nxy1=Nvec(iter);

    [Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,ldim,h,BCsPhi,Nmodes,"none",false) ;
    if iter ==1
        Phiref=Phi;
        Nxref=Nx;
        Nyref=Ny;
    
    else
        
        [Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om);

        %[Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om);
        
         Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly);
        %Phiref=Phi;
        
    end

    Nxy2=Nx*Ny;

    [Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,ldim,h,BCsPsi,Nmodes,"none",false) ;
     if iter ==1
        Psiref=Psi;
    else
        
        [Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly,Om2);

        %[Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly,Om2);

        Psi = eigensign(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly);

        % Psiref=Psi;
        % Nxref=Nx;
        % Nyref=Ny;
    end

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

                                %Hmat=Psik.*LPhipPhiq;
                                Hmat=Phip;
                                mdShape = reshape(Hmat,[(Ny+1),(Nx+1)]) ;
                                xax = (0:Nx)*h ;
                                yax = (0:Ny)*h ;
                                [X,Y] = meshgrid(xax,yax) ;


                                % subplot(3,3,ite2)
                                % pcolor(X,Y,3000*(mdShape));%,(abs(mdShape)));%,'FaceColor','texturemap') ;
                                % shading interp
                                if ite1 <17 
                                    ite1=ite1+1;
                                    ite2=ite2+1;
                                end
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
    
    H111(iter)=Hv(ntest,5,5);
    H121(iter)=Hv(ntest,6,5);
    H131(iter)=Hv(ntest,7,5);
    H141(iter)=Hv(ntest,8,5);
    H211(iter)=Hv(ntest,5,6);
    H221(iter)=Hv(ntest,6,6);
    H231(iter)=Hv(ntest,7,6);
    H241(iter)=Hv(ntest,8,6);
    H311(iter)=Hv(ntest,5,7);
    H321(iter)=Hv(ntest,6,7);
    H331(iter)=Hv(ntest,7,7);
    H341(iter)=Hv(ntest,8,7);
    H411(iter)=Hv(ntest,5,8);
    H421(iter)=Hv(ntest,6,8);
    H431(iter)=Hv(ntest,7,8);
    H441(iter)=Hv(ntest,8,8);

    sumH(1,iter)=sum(sum(sum(Hv(1:Nmodes,:,:))));
    
    % Hsp=Hv;
    % Hsp(abs(Hv)<1e-1)=0;
    % switch iter
    %     case 10
    %         Hsp10=Hsp;
    %     case 9
    %         Hsp9=Hsp;
    %     case 8
    %         Hsp8=Hsp;
    %     case 7
    %         Hsp7=Hsp;
    % 
    %     case 6
    %         Hsp6=Hsp;
    % 
    % 
    % end
    
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
% subplot(2,2,1)
semilogx(hvec,sumH(1,:),LineWidth=3,Marker="o")
xlabel('h')
ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^1_{i,j}'])
set (gca,'xdir','reverse')
set(gca,'FontSize',20)
% 
% subplot(2,2,2)
% semilogx(hvec,sumH(11,:),LineWidth=3,Marker="o")
% xlabel('h')
% ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^2_{i,j}'])
% set (gca,'xdir','reverse')
% set(gca,'FontSize',20)
% 
% subplot(2,2,3)
% semilogx(hvec,sumH(12,:),LineWidth=3,Marker="o")
% xlabel('h')
% ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^3_{i,j}'])
% set (gca,'xdir','reverse')
% set(gca,'FontSize',20)
% 
% subplot(2,2,4)
% semilogx(hvec,sumH(13,:),LineWidth=3,Marker="o")
% xlabel('h')
% ylabel(['Value of \Sigma_{i,j}^{' num2str(Nmodes) '} H^4_{i,j}'])
% set (gca,'xdir','reverse')
% set(gca,'FontSize',20)

%%
% Htest=sparse(squeeze(Hsp9(10,:,:)))
% figure
% spy(Htest)
% title('Non zero H^{10}_{i,j} coefficients for fully free plate')
% xlabel('i modes')
% ylabel('j modes')
% set(gca,'FontSize',20)
%%

% test=Hsp10-Hsp9;
% test(abs(test)<1e2)=0;
% 
% test=sparse(squeeze(test(1,:,:)));
% figure
% spy(test)
%% Save parameters

if ~exist("./param/", 'dir')
       mkdir("./param/")
end

save('./param/Test_NL_Fullclamp_7.mat','rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv');