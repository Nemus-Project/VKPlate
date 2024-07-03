clear all
close all
clc
%%%
%--- physical and elastic parameters

load("./Test_NL_Fullclamp_3.mat");  % loading parameters

%--- derived parameters (don't change here)
D       = E * Lz^3 / 12 / (1-nu^2) ;
k       = (1/44100)/2 ;             % Timestep
T       = 2;                        % Total time
Ts      = floor(T/k) ;              % Number of time grid points
t       =linspace(0,T,Ts);          % Time vector
x       =linspace(0,Lx,Nx+1);       % x space vector
y       =linspace(0,Ly,Ny+1);       % y space vector


%--- mode reconstruction
Fext=zeros(1,Ts);
c0=D/(rho*Lz);
T0=0.0005;
Thw=0.0005;
for n=1:Ts
    if abs(n*k-T0) <Thw
        Fext(1,n)=0.5*c0*(1+cos(pi*(n*k-T0)/Thw));
    end
end
figure
plot(t,Fext)




%colormap('parula') ;
for m = 1 : Nmodes
    mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
    %subplot(4,4,m)
    %mesh(3000*(mdShapes(:,:,m)),(abs(mdShapes(:,:,m))),'FaceColor','texturemap') ;  %axis equal; axis tight ; view(2) ;
    pext(:,m)=0e-1*Fext*mdShapes(floor(end/4),floor(end/4),m);

end
pext=pext';

%%
%chi=0.1*sort(1+1*rand(Nmodes,1));
Id=speye(Nmodes);                               % Identity matrix

Dw1=diag((Om(1:Nmodes)));                       % Diagonal matrix of Omega

Dw2=diag((Om(1:Nmodes)).^2);                    % Diagonal matrix of Omega^2

Dw=diag((cos(k*Om(1:Nmodes))));                 % Diagonal for the exact lossless integrator

DV0=2*k^(-2)*(Id-diag((cos(k*Om(1:Nmodes)))));  % 

Ddiag=2*sparse(Dw);                             % 2*Sparse of Dw

D0sv=2*Id-(k^2)*Dw2;                            % 

%%
qm=zeros(Nmodes,1);                             % 
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
etam=zeros(Nmodes,1);
qh=zeros(Nmodes,Nmodes,Nmodes);
%qm=flipud(sort(1e-4*rand(Nmodes,1)));      %random excitation
xi=1e-1*Id;
DSM=(0.5*k*xi*Dw1+Id)^-1;
Dqmsv=0.5*k*xi*Dw1-Id;
qm(1)=5e-3;
q0=qm;

qmsv=qm;
q0sv=qmsv;

Q(1,:)=q0;
Qsv(1,:)=q0sv;
enerr=1/(rho*Lz);

%% Initialisation


for ki = 1:Nmodes
    for ji =1:Nmodes
        for ii =1:Nmodes
            eta(ki)=eta(ki)-((E*Lz)/(2*zetafourth(ki)))*Hv(ii,ji,ki)*q0(ii)*q0(ji);

        end
    end
end
etasv=eta;

ETA(1,:)=eta;
ETASV(1,:)=etasv;

U=enerr*(1/(2*E*Lz))*zetafourth(1:Nmodes)'*0.5*(eta.^2+eta.^2);   %to change with actual energy (meaning etan-1)
%U=enerr*(1/(2*E*Lz))*zetafourth(1)*0.5*(eta.^2+eta.^2);

ep=1e-10;
psim1=sqrt(2*U+ep);

%%


%w=zeros(Ny+1,Nx+1);
tic
for n = 2 : Ts

    En(n)=0.5*sum(((q0-qm)/k).^2);

    V0(n)=0.5*sum(DV0*q0.*qm);

    eta(:)=0;
    etasv(:)=0;
    for li = 1:Nmodes
        for mi =1:Nmodes
            for ni =1:Nmodes
                eta(li)=eta(li)-((E*Lz)/(2*zetafourth(li)))*Hv(mi,ni,li)*q0(mi)*q0(ni);
                etasv(li)=etasv(li)-((E*Lz)/(2*zetafourth(li)))*Hv(mi,ni,li)*q0sv(mi)*q0sv(ni);

            end
        end
    end

    ETA(n,:)=eta;

    nlsv(:)=0;

    for ki = 1:Nmodes
        for ji =1:Nmodes
            for ii =1:Nmodes
                nlsv(ki)=nlsv(ki)+Hv(ki,ji,ii)*q0sv(ii)*etasv(ji);



            end
        end
    end
    % for si = 1 : Nmodes
    %     for ki = 1 : Nmodes
    %     for li = 1:Nmodes
    %     for mi =1:Nmodes
    %         for ni =1:Nmodes
    %
    %             nlsv(si)=nlsv(si)-Hv(li,mi,ni)*Hv(ki,si,li)*q0sv(ki)*q0sv(mi)*q0sv(ni)/(2*zetafourth(li));
    %
    %
    %         end
    %     end
    %     end
    %     end
    % end



    U=enerr*(1/(2*E*Lz))*zetafourth(1:Nmodes)'*(eta.^2);
    %U=enerr*(1/(2*E*Lz))*zetafourth(1)*(eta.^2);


    dUdq(:)=0;
    for ji = 1:Nmodes
        for ki =1:Nmodes
            for ii =1:Nmodes
                dUdq(ji)=dUdq(ji)-eta(ki)*(q0(ii)*Hv(ii,ji,ki));

            end
        end
    end
    %dUdq=-2*dUdq;
    g=enerr*(1/sqrt(2*U+ep))*dUdq;
    G(n)=g(1);

    Q(n,:)=q0;
    Qsv(n,:)=q0sv;


    Dqnp=0.25*(k^2)*g*g'+Id;
    Dqnm=0.25*(k^2)*g*g'-Id;
    %Dexp=inv(Dqnp);
    Dexp=Id-(0.25*(k^2)*g*g')/(1+0.25*k^2*g'*g);

    Dqnm2=0.25*(k^2)*g*g'+0.5*xi*Dw1*k-Id;
    %Dqnp2=0.25*(k^2)*g*g'+0.5*xi*Dw1*k+Id;

    %Dexp2=inv(Dqnp2);
    Dexp2=DSM-(DSM*0.25*(k^2)*g*g'*DSM)/(1+0.25*k^2*g'*DSM*g);

    qs=Dexp2*(D0sv*q0+Dqnm2*qm-g*(k^2)*psim1);% +k^2*pext(:,n));

    %qs=Dexp*(Ddiag*q0+Dqnm*qm-g*(k^2)*psim1);%no damping


    % qssv=D0sv*q0sv-qmsv+((k^2)/(rho*Lz))*nlsv;%no damping
    qssv=DSM*(D0sv*q0sv+Dqmsv*qmsv+((k^2)/(rho*Lz))*nlsv);
    %qssv=D0sv*q0sv-qmsv+((k^2)*E/(rho))*nlsv;




    psi0=psim1+0.5*g'*(qs-qm);




    % w(:,:)=0;
    % for m = 1 : Nmodes
    %     w=w+q0(m)*mdShapes(:,:,m);
    % end
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
scatter(t,(EN-EN(150))/EN(150),'.')
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
figure
plot(Q-Qsv,LineWidth=3)
title('w_k error')
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

%%
figure
plot(PSI0,LineWidth=3)
title('Psi')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
