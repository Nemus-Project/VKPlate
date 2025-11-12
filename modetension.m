%%---- test VK coefficient calc
% What is This?
% and what is it meant to do?
% who wrote it?
clear all
close all
clc
addpath '/Users/alexis/Documents/MATLAB/VKPlate/private/magpie'
colorMap = load('Ha3ColorMap.mat').custom_map;


%% ------------------------------------------------------------------------
% custom params
% rho     = 7860 ;
% E       = 2e+11 ;
%nu      = 0.3 ;
E       = 4e+9 ;
rho     = 900 ; 
nu      = 0.4 ;
Lz      = 1e-5 ;
Lx      = 4.3e-2 ;
Ly      = Lx ;
T       = 0 ;
Nmodes  =10;
Npsi=50;
npts=8;
%hvec=[0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003,0.002];
%Nvec=floor(logspace(1.5,2.5,npts));
Nvec=[31 38 44 51 61 71 85 100];
tmagpie2=zeros(npts,1);
tvk=zeros(npts,1);
%Nvec=floor(logspace(2,3.2,10));
%Nvec=[1000:100:2000];
%Hconv=zeros(Nmodes,Nmodes,Nmodes,length(Nvec));
compnorm=zeros(1,length(Nvec));

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
nxtestvec=ones(npts,1);
ntest=1;
ktest=ntest;
ptest=1;
qtest=1;
ltest=1;
mtest=1;
ite1=2;
ite2=1;
figure
Gamcon=zeros(Nmodes,npts);%convergence matrix
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

    [Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,T,ldim,h,BCsPhi,Nmodes,"none",true) ;%calls magpie
     Nxphi(iter)=Nx;

  
   for nQ = 1 : Nmodes % normalize the basis against the first value
        if  nQ == 1
        Phitemp   = Phi(:,nQ) ;
        Phinorm   = trapzIntcalc(Phitemp.*Phitemp,h,Nx,Ny) ;
        Phitemp   = Phitemp / sqrt(Phinorm) ;
        Phi(:,nQ) = Phitemp ;
        else
        Phitemp   = Phi(:,nQ) ;
        Phitemp   = Phitemp / sqrt(Phinorm) ;
        Phi(:,nQ) = Phitemp ;
        end
    end
    compnorm(iter)=Phinorm;

    % [Phiort,OrtN]=GramSchmidt(Phi); %orthogonalise the basis of phi
    % 
    % for nQ = 1 : Nmodes % normalize the basis
    %     Phitemp   = Phiort(:,nQ) ;
    %     Phinorm   = trapzIntcalc(Phitemp.*Phitemp,h,Nx,Ny) ;
    %     Phitemp   = Phitemp / sqrt(Phinorm) ;
    %     Phi(:,nQ) = Phitemp ;
    % end
    % clear Phitemp
    % 
    % 
    % 
    % if iter ==1
    %     Phiref=Phi;
    %     Nxref=Nx;
    %     Nyref=Ny;
    % 
    % else
    % 
    %     [Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om);
    % 
    %     %[Phiort,OrtN]=GramSchmidt(Phi); %orthogonalise the basis of phi,
    %     %not necessarry, just the normalize is useful
    % 
    %     for nQ = 1 : Nmodes % normalize the basis
    %         Phitemp   = Phi(:,nQ) ;
    %         Phinorm   = trapzIntcalc(Phitemp.*Phitemp,h,Nx,Ny) ;
    %         Phitemp   = Phitemp / sqrt(Phinorm) ;
    %         Phi(:,nQ) = Phitemp ;
    %     end
    %     clear Phitemp
    % 
    % 
    %     %[Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om);
    % 
    %     Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,T,Nmodes,Lx,Ly);
    %     %Phiref=Phi;
    % 
    % end

    Nxy2=Nx*Ny;

    [Om2,Psi,Nxp,Nyp,~,zetafourth] = magpie(rho,E,nu,T,ldim,h,BCsPsi,Npsi,"none",true) ;
     Nxpsi(iter)=Nxp;


    for nQ = 1 : Nmodes % normalize the basis against the first value
        if  nQ == 1
        Psitemp   = Psi(:,nQ) ;
        Psinorm   = trapzIntcalc(Psitemp.*Psitemp,h,Nx,Ny) ;
        Psitemp   = Psitemp / sqrt(Psinorm) ;
        Psi(:,nQ) = Psitemp ;
        else
        Psitemp   = Psi(:,nQ) ;
        Psitemp   = Psitemp / sqrt(Phinorm) ;
        Psi(:,nQ) = Psitemp ;
        end
    end
    % [Psiort,OrtN]=GramSchmidt(Psi); %orthogonalise the basis of phi
    % 
    % for nQ = 1 : Nmodes % normalize the basis
    %     Psitemp   = Psiort(:,nQ) ;
    %     Psinorm   = trapzIntcalc(Psitemp.*Psitemp,h,Nx,Ny) ;
    %     Psitemp   = Psitemp / sqrt(Psinorm) ;
    %     Psi(:,nQ) = Psitemp ;
    % end
    % clear Phitemp
    % 
    % 
    % if iter ==1
    %     Psiref=Psi;
    % else
    % 
    %     [Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Npsi,Lx,Ly,Om2);
    % 
    % 
    %     %[Psiort,OrtN]=GramSchmidt(Psi); %orthogonalise the basis of phi
    %     %not necessarry, just the normalize is useful
    % 
    %     for nQ = 1 : Nmodes % normalize the basis
    %         Psitemp   = Psi(:,nQ) ;
    %         Psinorm   = trapzIntcalc(Psitemp.*Psitemp,h,Nx,Ny) ;
    %         Psitemp   = Psitemp / sqrt(Psinorm) ;
    %         Psi(:,nQ) = Psitemp ;
    %     end
    %     clear Phitemp
    % 
    %     %[Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly,Om2);
    % 
    %     Psi = eigensign(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Npsi,Lx,Ly);
    % 
    %     % Psiref=Psi;
    %     % Nxref=Nx;
    %     % Nyref=Ny;
    % end

    disp(iter)
    %assert((Nxy1-Nxy2 ==0),'Nxy pure man facked')

    zeta = (zetafourth).^(1/4) ;

    Hv = zeros(Npsi,Ntensor,Ntensor) ;
    %Ev = zeros(Ntensor,Ntensor,Ntensor) ;

    Dxx = DxxBuild(Nx,Ny,h) ;
    Dyy = DyyBuild(Nx,Ny,h) ;
    Dxy = DxyBuild(Nx,Ny,h) ;

    tic;
    for k = 1 : Npsi
        %Phik = Phi(:,k) ;
        Psik = Psi(:,k) ;
        %Phiknorm   = trapzIntcalc(Phik.*Phik,h,Nx,Ny);
        %Psiknorm   = trapzIntcalc(Psik.*Psik,h,Nx,Ny);
        for p = 1 : Ntensor
            Phip = Phi(:,p);
            Phipnorm   = trapzIntcalc(Phip.*Phip,h,Nx,Ny) ;
            for q = p : Ntensor
                Phiq = Phi(:,q) ; Psiq = Psi(:,q);

                %Phiqnorm   = trapzIntcalc(Phiq.*Phiq,h,Nx,Ny) ;
                %Psiqnorm   = trapzIntcalc(Psiq.*Psiq,h,Nx,Ny) ;

                LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;

                %LPhipPsiq = vkOperator(Phip,Psiq,Dxy,Dxx,Dyy) ;

                Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);%/sqrt(Psiknorm*Phipnorm*Phiqnorm) ;
                %Hv(p,k,q) = Hv(k,p,q);
                %Hv(k,q,p) =  Hv(k,p,q);
                Hv(k,q,p) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);%/sqrt(Psiknorm*Phipnorm*Phiqnorm);
                %Ev(k,p,q) = trapzIntcalc(Phik.*LPhipPsiq,h,Nx,Ny);%/sqrt(Phiknorm*Phipnorm*Psiqnorm) ;

                if iter == ite1

                    if k==ktest

                        if p==ptest

                            if q==qtest


                                %Hmat=Psik.*LPhipPhiq;
                                % Hmat=Phip;
                                % mdShape = reshape(Hmat,[(Ny+1),(Nx+1)]) ;
                                xax = (0:Nx)*h ;
                                yax = (0:Ny)*h ;
                                [X,Y] = meshgrid(xax,yax) ;

                                %
                                % subplot(3,3,ite2)
                                % %figure
                                % pcolor(X,Y,3000*(mdShape));%,(abs(mdShape)));%,'FaceColor','texturemap') ;
                                % daspect([1 1 Lx/Ly])
                                % set(gca,'xtick',[])
                                % set(gca,'ytick',[])
                                % title([ num2str(Ny+1) '\times' num2str(Nx+1)])
                                % set(gca,'Fontsize',20)
                                % colormap(colorMap)
                                % shading interp
                                % if ite1 <17
                                %     ite1=ite1+1;
                                %     ite2=ite2+1;
                                % end
                            end
                        end
                    end


                end
            end
        end
    end
    tvk=toc;
    %  mdShape1 = reshape(Hmat(2,2,1,:),[(Ny+1),(Nx+1)]) ;
    % % clear Hmat
    %               xax = (0:Nx)*h ;
    %              yax = (0:Ny)*h ;
    %              [X,Y] = meshgrid(xax,yax) ;
    % %
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

    sumH(1,iter)=sum(sum(sum(abs(Hv(1:Nmodes,:,:)))));

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
    % if Nx*h-Nxref*hvec(1) ~= 0
    %     disp("x is fudged")
    %     nxtestvec(iter)=0;
    %     %fix(Lx/h)
    %     % Lx=Lx+h/2;
    %     % Ly=Ly+h/2;
    % end
    % if Ny*h-Nyref*hvec(1) ~= 0
    %     disp("y is fudged")
    %     nxtestvec(iter)=0;
    %     %fix(Lx/h)
    %     % Lx=Lx+h/2;
    %     % Ly=Ly+h/2;
    % end
    clear Hv


    %assert((Nx*h-Nxref*hvec(1)==0),'ALERT ROUNDING ERROR ALERT')
    %toc
    %nvec(iter)=Nx*Ny;

    Ev=permute(Hsp,[3 1 2]);fi


    filcoeff=0;

    Npsivec=[1:Npsi];
    p=5;

    for l = 1 : Npsi

        filcoeff=filcoeff+ real(Hsp(l,p,p)*Ev(p,l,p))/(2*zeta(l)^4);
        Gamcon(l,iter) = filcoeff;

    end
    filcoeff=0;
%Tests
if iter==1
    OmPsi1=Om2;
else
    OmPsi2=Om2;
end

end
plotnum=find(nxtestvec);
Hv=Hsp;
%%

[Xconv,Yconv] = meshgrid(hvec(1:npts)',Npsivec(1:Npsi));
figure
surf(Xconv',Yconv',real(Gamcon(1:Npsi,1:npts)'))
shading interp
xlabel("Grid spacing (m)")
ylabel("Number of in plane modes N_\Psi")
zlabel("\Gamma_{ppp}^p")
hold on
hplot=npts-1;
Nvecplot=Npsi-10;
[Xplt1,Yplt1] = meshgrid(hvec(hplot)*ones(length(hvec))',Npsivec(1:Npsi));
plot3(Xplt1',Yplt1',squeeze(real(Gamcon(1:Npsi,hplot))'),LineWidth=4,Color='k')

[Xplt2,Yplt2] = meshgrid(hvec(1:npts)',Npsivec(Nvecplot)*ones(Npsi));
plot3(Xplt2',Yplt2',squeeze(real(Gamcon(Nvecplot,1:npts))'),LineWidth=4)
%p.LineWidth=4;
%set(gcf,"LineWidth",6)
set(gca,"Fontsize",26)
%colormap gray
%%
figure
semilogx(hvec(1,plotnum),squeeze(real(Gamcon(Nvecplot,plotnum))),LineWidth=3,Marker="o")
set (gca,'xdir','reverse')
set(gca,'FontSize',20)
ylabel("\Gamma_{ppp}^p")
xlabel("Grid spacing (m)")
figure
plot(Npsivec(1:Npsi),squeeze(real(Gamcon(1:Npsi,hplot))),LineWidth=3,Color='k')
%set (gca,'xdir','reverse')
set(gca,'FontSize',20)
ylabel("\Gamma_{ppp}^p")
xlabel("Number of in plane modes N_\Psi")

%%



%%


%%
Htest=sparse(squeeze(Hsp(1,:,:)));
figure
spy(Htest)
title('Non zero H^{1}_{i,j} coefficients for fully free plate')
xlabel('i modes')
ylabel('j modes')
set(gca,'FontSize',20)
%% Orthomap
Orthog=zeros(Nmodes);
for i=1:Nmodes
    for j=1:Nmodes
        tem=Phi(:,i)'*Phi(:,j);
        if abs(tem)<1e-5
            tem=0;
        else
            tem=1;
        end
        Orthog(i,j)=tem;
    end
end
figure
imagesc(Orthog)
title("Phi")
colormap(flipud(gray))

%% Orth
[Phiort,OrtN]=GramSchmidt(Phi);
Orthog=zeros(Nmodes);
for i=1:Nmodes
    for j=1:Nmodes
        tem=Phiort(:,i)'*Phiort(:,j);
        if abs(tem)<1e-5
            tem=0;
        else
            tem=1;
        end
        Orthog(i,j)=tem;
    end
end

figure
imagesc(Orthog)
title("post GS")
colorbar
colormap(flipud(gray))
%%

%%

%%
mus=1;
xax = (0:Nx)*h ;
yax = (0:Ny)*h ;
[X,Y] = meshgrid(xax,yax) ;
mdShape = reshape(Phi(:,mus),[(Ny+1),(Nx+1)]);
ma=max(max(abs(mdShape)));
figure%('Renderer', 'Painters')
surf(X,Y,-(mdShape)./ma);
%alpha 0.9
%daspect([1 1 Lx/Ly])
%colorMap = load('Ha3ColorMap.mat').custom_map;
%colormap(colorMap)
set(gca,'XTick',[])
set(gca,'YTick',[])
%colormap(randmap)
%caxis([-ma ma])
caxis([-1 1])
shading interp
%colorbar
%% Test Tension modeshape
close all
Phimem=zeros(size(mdShape));
Phimem=sin(pi*xax/Lx)'*sin(pi*yax/Ly);

Phimemlin = reshape(Phimem,[(Ny+1)*(Nx+1),1]);
normem=trapzIntcalc(Phimemlin.*Phimemlin,h,Nx,Ny)
normpla=trapzIntcalc(Phi(:,mus).*Phi(:,mus),h,Nx,Ny)

Phimem2=Phimem/sqrt(normem);

figure
surf(Phimem2)
figure
surf(mdShape(:,:))

Phierr=100*abs(Phimem2(:,:)-abs(mdShape(:,:)))/max(max(Phimem2(:,:)));

figure('Renderer', 'Painters')
pcolor(Phierr)
daspect([1 1 Lx/Ly])
set(gca,'XTick',[])
set(gca,'YTick',[])
caxis([0 1])
colorbar
shading interp
set(gca,'FontSize',24)
