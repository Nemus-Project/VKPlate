clear all
close all
clc 
%%
%--- physical and elastic parameters
%h       = sqrt(Lx*Ly)*0.01;
load("./Test_NL_Fullclamp_3.mat");

%--- parameters 
k       = 1/44100 ;% Timestep
T       = 0.1;
Ts      = floor(T/k) ;% Number of time grid points
t       =linspace(0,T,Ts);
x       =linspace(0,Lx,Nx+1);
y       =linspace(0,Ly,Ny+1);

%--- mode reconstruction
Nmodes=1;
%%
Id=1;

Dw2=Om(1)^2;

Dw=(cos(k*Om(1)));

DV0=2*k^(-2)*(Id-cos(k*Om(1)));

Ddiag=2*Dw;

%%
qm=zeros(Nmodes,1);

q0=zeros(Nmodes,1);

qmsv=zeros(Nmodes,1);
q0sv=zeros(Nmodes,1);


dUdq=zeros(Nmodes,1);
g=zeros(Nmodes,1);

En=zeros(1,Ts);

V0=zeros(1,Ts);

PSI0=zeros(1,Ts);

Q=zeros(Ts,Nmodes);

QSV=zeros(Ts,Nmodes);

ETA=zeros(Ts,Nmodes);

ETASV=zeros(Ts,Nmodes);

eta=zeros(Nmodes,1);

qh=zeros(Nmodes,Nmodes,Nmodes);
%% Initial displacement
qm(1)=1e-3;
%%
q0=qm;

q0sv=q0;

qmsv=qm;

Q(1,:)=q0;

QSV(1,:)=q0sv;

%% Initialisation



eta=eta-((E*Lz)/(2*zetafourth(1)))*Hv(1,1,1)*q0*q0;

etasv=eta;
   
ETA(1,:)=eta;
ETASV(1,:)=etasv;

U=(1/(2*E*Lz))*zetafourth(1)*0.5*(eta.^2+eta.^2);   %to change with actual energy

ep=1e-14;
psim1=sqrt(2*U+ep);

%%



tic
for n = 2 : Ts

    En(n)=0.5*sum(((q0-qm)/k).^2);

    V0(n)=0.5*sum(DV0*q0.*qm);
    
    eta=-((E*Lz)/(2*zetafourth(1)))*Hv(1,1,1)*q0(1)*q0(1);

    etasv=-((E*Lz)/(2*zetafourth(1)))*Hv(1,1,1)*q0sv(1)*q0sv(1);
                
    ETA(n,:)=eta;

    ETASV(n,:)=etasv;

    U=(1/(2*E*Lz))*zetafourth(1)*(eta.^2);
    
    dUdq=-eta*(q0*Hv(1,1,1));
           
    g=(1/sqrt(2*U+ep))*dUdq;
    
    G(n)=g;
    
    Q(n,:)=q0;

    QSV(n,:)=q0sv;
    

    %qs=(1/(1+(g^2*k^2)/4))*((2*cos(Om(1)*k))*q0+(((g^2*k^2)/4)-1)*qm-k^2*g*psim1);%Exact integrator
    
    qs=(1/(1+(g^2*k^2)/4)) * ((2-Om(1)^2*k^2)*q0+(((g^2*k^2)/4)-1)*qm-k^2*g*psim1); %SAV but not exact

    qssv=(2-k^2*Om(1)^2)*q0sv - qmsv + ((Hv(1,1,1)*k^2)/(rho*Lz))*q0sv*etasv ;%Stormer Verlet integrator

    psi0=psim1+0.5*g*(qs-qm);

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
scatter(t,(EN-EN(150))./EN,".")
xlim([0.02 T])

%% 
figure
plot(t,ETA(:,:),t,ETASV(:,:),LineWidth=3)
title('\eta_k')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
figure
plot(t,Q(:,:),t,QSV(:,:),LineWidth=3)
title('w_k SAV & Stormer Verlet')
xlabel('Times')
ylabel('Amplitude')
set(gca,'FontSize',20)


%%
err=(Q(:,:)-QSV(:,:));
figure
plot(err,LineWidth=3)
title('Absolute error')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)


