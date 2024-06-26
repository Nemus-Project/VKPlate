clear all
close all
clc
%%
%--- physical and elastic parameters

load("./Test_NL_Fullclamp_4.mat");

%--- normalisation of the coupling coefficient (useless)
%Hv=abs(Hv);
%mHv=max(max(max(Hv)));
%Hv=10*abs(Hv)/mHv;
%Hv(abs(Hv)<1e-5)=0;

%--- derived parameters (don't change here)
D       = E * Lz^3 / 12 / (1-nu^2) ;
k       = (1/44100)/4 ;% Timestep
T       = 0.1;
Ts      = floor(T/k) ;% Number of time grid points
t       =linspace(0,T,Ts);
x       =linspace(0,Lx,Nx+1);
y       =linspace(0,Ly,Ny+1);

%--- mode reconstruction

for m = 1 : Nmodes
    mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;

end

%%
chi=0.1*sort(1+1*rand(Nmodes,1));
Id=speye(Nmodes);
Dw2=diag((Om).^2);
Dw=diag((cos(k*Om)));
DV0=2*k^(-2)*(Id-diag((cos(k*Om))));
Ddiag=2*sparse(Dw);
Dssv=2*Id-(k^2)*Dw2;

%%
qm=zeros(Nmodes,1);
q0=zeros(Nmodes,1);
nlsv=zeros(Nmodes,1);

qmsv=zeros(Nmodes,1);
q0sv=zeros(Nmodes,1);
qssv=zeros(Nmodes,1);

dUdq=zeros(Nmodes,1);
g=zeros(Nmodes,1);
En=zeros(1,Ts);
V0=zeros(1,Ts);
PSI0=zeros(1,Ts);
Q=zeros(Ts,Nmodes);
Qsv=zeros(Ts,Nmodes);
ETA=zeros(Ts,Nmodes);
%W=zeros(Ny+1,Nx+1,Ts);
eta=zeros(Nmodes,1);
etasv=zeros(Nmodes,1);
etam=zeros(Nmodes,1);
qh=zeros(Nmodes,Nmodes,Nmodes);
%qm=flipud(sort(1*rand(Nmodes,1)));
qm(1)=0.001;
q0=qm;

qmsv=qm;
q0sv=qmsv;

%% Initialisation


for ki = 1:Nmodes
    for ji =1:Nmodes
        for ii =1:Nmodes
            eta(ki)=eta(ki)-((E*Lz)/(2*zetafourth(ki)))*Hv(ii,ji,ki)*q0(ii)*q0(ji);

        end
    end
end


ETA(1,:)=eta;
U=(1/(2*E*Lz))*zetafourth'*(eta.^2+eta.^2);   %to change with actual energy

ep=1e-10;
psim1=sqrt(2*U+ep);

%%


w=zeros(Ny+1,Nx+1);
tic
for n = 2 : Ts

    En(n)=0.5*sum(((q0-qm)/k).^2);
    V0(n)=0.5*sum(DV0*q0.*qm);




    etam(:)=0;
    eta(:)=0;
    etasv(:)=0;
    for ki = 1:Nmodes
        for ji =1:Nmodes
            for ii =1:Nmodes
                etam(ki)=etam(ki)-((E*Lz)/(2*zetafourth(ki)))*Hv(ii,ji,ki)*qm(ii)*qm(ji);
                eta(ki)=eta(ki)-((E*Lz)/(2*zetafourth(ki)))*Hv(ii,ji,ki)*q0(ii)*q0(ji);
                %etasv(ki)=etasv(ki)-((E*Lz)/(2*zetafourth(ki)))*Hv(ii,ji,ki)*q0sv(ii)*q0sv(ji);

            end
        end
    end


    nlsv(:)=0;
    for ki = 1:Nmodes
        for ji =1:Nmodes
            for ii =1:Nmodes
                nlsv(ki)=nlsv(ki)+Hv(ki,ji,ii)*q0sv(ii)*eta(ji);

            end
        end
    end

    ETA(n,:)=eta;
    U=(1/(2*E*Lz))*zetafourth'*(eta.^2);

    dUdq(:)=0;
    for ji = 1:Nmodes
        for ki =1:Nmodes
            for ii =1:Nmodes
                dUdq(ji)=dUdq(ji)-0.5*eta(ki)*(q0(ii)*Hv(ii,ji,ki));

            end
        end
    end
    %dUdq=-2*dUdq;
    g=(1/sqrt(2*U+ep))*dUdq;
    G(n)=g(1);

    Q(n,:)=q0;
    Qsv(n,:)=q0sv;


    Dqnp=0.25*(k^2)*g*g'+Id;
    Dqnm=0.25*(k^2)*g*g'-Id;
    Dexp=inv(Dqnp);
    %Dexp=
    qs=Dexp*(Ddiag*q0+Dqnm*qm-g*(k^2)*psim1);

    qssv=Dssv*q0-qm+((k^2)/(rho*Lz))*nlsv;



    psi0=psim1+0.5*g'*(qs-qm);




    w(:,:)=0;
    for m = 1 : Nmodes
        w=w+q0(m)*mdShapes(:,:,m);
    end
    %W(:,:,n)=w;
    PSI0(n)=psim1;
    qm=q0;
    q0=qs;
    psim1=psi0;




    qmsv=q0sv;
    q0sv=qssv;


end
toc











%%
figure
plot(t,En-En(3),t,V0-V0(3),t,((PSI0.^2)-(PSI0(3)^2))/2)
figure
EN=En+V0+(PSI0.^2)/2;
plot(t,(EN-EN(150))/EN(150))
xlim([0.02 T])
%%
% Fs=44100;
% solo=squeeze(W(floor(Nx/3),floor(Ny/4),:));
% fftest=fft(solo);
% P2 = (fftest/Ts);
% P1 = P2(1:Ts/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(Ts/2))/Ts;
% fftsig=P1;
%
% figure
% plot(f,P1)
% hold on
% xline(Om)
% xlim([0 Om(end)])

%% Anim

% az=00;
% el=-70;
% mmoutu=max(max(max(W(:,:,:))));
% figure('Position', [1000 1000 600 600],'Color','w')
% surf(x,y,W(:,:,1),"LineWidth",3)
% %view(az,el)
% %light("Style","local","Position",[-10 -10 0]);
% %lighting gouraud;material dull ;
% set(gca,'Color','w');grid off;
% colormap(HA2);shading interp;
% zlim([-mmoutu 1*mmoutu])
% caxis([-mmoutu mmoutu])
% xticks([]);yticks([]);zticks([]);
% axis off
% drawnow
% gif('modeplatemicheleup.gif')
% for n=11:10:500
% surf(x,y,W(:,:,n),"LineWidth",3)
% %light("Style","local","Position",[-10 -10 0]);
% %lighting gouraud
% set(gca,'Color','w')
% %view(az,el)
% caxis([-mmoutu mmoutu])
% material dull
% grid off;xticks([]);yticks([]);zticks([]);
% axis off
% %colormap(HA2);
% shading interp
% zlim([-mmoutu 1*mmoutu])
% drawnow
% gif
% pause(0.05)

%end
%%
figure
plot(ETA(:,:),LineWidth=3)
title('\eta_k')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
figure
plot(Q(:,:),LineWidth=3)
title('w_k SAV')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
figure
plot(Qsv(:,:),LineWidth=3)
title('w_k Stormer-Verlet')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)


%%
% 
% figure
% for anim=2:5:length(t)
%     surf(W(:,:,anim))
%     shading interp
%     zlim([-1000 1000])
%     pause(0.01)
% 
% end

