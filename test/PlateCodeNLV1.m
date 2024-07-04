clear all
close all
clc 
%%
%--- physical and elastic parameters
%h       = sqrt(Lx*Ly)*0.01;
load("./param/modesteswb.mat");
HA2=load("Ha2ColorMap.mat").Ha2ColorMap;


%--- derived parameters (don't change here)
D       = E * Lz^3 / 12 / (1-nu^2) ;
k       = 1/44100 ;% Timestep
T       = 1.5;
Ts      = floor(T/k) ;% Number of time grid points
t       =linspace(0,T,Ts);
x       =linspace(0,Lx,Nx+1);
y       =linspace(0,Ly,Ny+1);

%--- mode reconstruction

%%
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
    mdShapes(:,:,m) = reshape(V(:,m),[(Ny+1),(Nx+1)]) ;
   %subplot(4,4,m)
   %mesh(3000*(mdShapes(:,:,m)),(abs(mdShapes(:,:,m))),'FaceColor','texturemap') ;  %axis equal; axis tight ; view(2) ;
   pext(:,m)=100*rho*h^(-3)*Fext*mdShapes(floor(end/4),floor(end/4),m);
   
end
pext=pext';
%pext(:,:)=0;
%%
chi=0.1*sort(1+1*rand(Nmodes,1));
Id=speye(Nmodes);
Dw2=diag((2*pi*freqs).^2);
Dw=diag((cos(k*2*pi*freqs)));
DV0=2*k^(-2)*(Id-diag((cos(k*2*pi*freqs))));
Ddiag=2*sparse(Dw);
%Ddiag=2*Id-k^2*Dw2;
Ddamp=2*diag(2*pi*freqs.*chi);
Dimpl=Id+k*Ddamp/2;
Dexp=sparse(inv(Dimpl));
%Dimpl=
Ddamp2=Id-k*Ddamp/2;

%%
qm=zeros(Nmodes,1);
q0=zeros(Nmodes,1);
En=zeros(1,Ts);
V0=zeros(1,Ts);
Q=zeros(Ts,Nmodes);
W=zeros(Ny+1,Nx+1,Ts);


% qm=flipud(sort(1*rand(Nmodes,1)));
% q0=qm;
%qm([1,5])=1;
%q0=qm;
w=zeros(Ny+1,Nx+1);

tic
for n = 1 : Ts

    En(n)=0.5*sum(((q0-qm)/k).^2);
    V0(n)=0.5*sum(DV0*q0.*qm);
    %qs=Ddiag*q0-qm-k*Ddamp*(q0-qm)+k^2*pext(:,n);
    %qs=Dimpl\(Ddiag*q0-Ddamp2*qm+k^2*pext(:,n));
    qs=Dexp*(Ddiag*q0-Ddamp2*qm+k^2*pext(:,n));
    qm=q0;
    q0=qs;
    
    Q(n,:)=q0;
    w(:,:)=0;
    for m = 1 : Nmodes
        w=w+q0(m)*mdShapes(:,:,m);
    end
    W(:,:,n)=w;
end
toc











%%
figure
plot(t,En-En(3))
hold on
plot(t,V0-V0(3))
figure
E=En+V0;
plot(t,(E-E(150))/E(150))
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
% xline(freqs)
% xlim([0 freqs(end)])

%%
az=00;
el=-70;
mmoutu=max(max(max(W(:,:,:))));
figure('Position', [1000 1000 600 600],'Color','w')
surf(x,y,W(:,:,1),"LineWidth",3)
%view(az,el)
%light("Style","local","Position",[-10 -10 0]);
%lighting gouraud;material dull ;
set(gca,'Color','w');grid off;
colormap(HA2);shading interp;
zlim([-mmoutu 1*mmoutu])
caxis([-mmoutu mmoutu])
xticks([]);yticks([]);zticks([]);
axis off
drawnow
gif('modeplatemicheleup.gif') 
for n=11:10:500
surf(x,y,W(:,:,n),"LineWidth",3)
%light("Style","local","Position",[-10 -10 0]);
%lighting gouraud
set(gca,'Color','w')
%view(az,el)
caxis([-mmoutu mmoutu])
material dull 
grid off;xticks([]);yticks([]);zticks([]);
axis off
%colormap(HA2);
shading interp
zlim([-mmoutu 1*mmoutu])
drawnow
gif
pause(0.05)

end

