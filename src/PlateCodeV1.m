clear all
close all
clc 
%%
%--- physical and elastic parameters
load("./param/modestest.mat");


%--- derived parameters (don't change here)
D       = E * Lz^3 / 12 / (1-nu^2) ;
k       = 1/44100 ;% Timestep
T       = 0.1;
Ts      = floor(T/k) ;% Number of time grid points
t       =linspace(0,T,Ts);
x       =linspace(0,Lx,Nx+1);
y       =linspace(0,Ly,Ny+1);
%--- mode reconstruction

%colormap('parula') ;
for m = 1 : Nmodes
    mdShapes(:,:,m) = reshape(V(:,m),[(Ny+1),(Nx+1)]) ;
   %subplot(4,4,m)
   %mesh(3000*(mdShapes(:,:,m)),(abs(mdShapes(:,:,m))),'FaceColor','texturemap') ;  %axis equal; axis tight ; view(2) ;
end

Id=speye(Nmodes);
Dw2=diag((2*pi*freqs).^2);
Dw=diag((cos(k*2*pi*freqs)));
DV0=2*k^(-2)*(Id-diag((cos(k*2*pi*freqs))));
Ddiag=2*sparse(Dw);



qm=zeros(Nmodes,1);
q0=zeros(Nmodes,1);
En=zeros(1,Ts);
V0=zeros(1,Ts);
Q=zeros(Ts,Nmodes);
W=zeros(Ny+1,Nx+1,Ts);
qm=sort(1*rand(Nmodes,1));
q0=qm;
w=zeros(Ny+1,Nx+1);

%psim1=u0*omeg0;
%Hsav(n-1)=(1/2)*((u0-um1)/k)^2+(1/2)*(psim1)^2;
%Hsav(n-1)=(1/(2*k^2))*(sum((u0(:)-um(:)).^2))+(1/2)*(psim1)^2;

%V0=(c^2/h^2)*(uh(1)^2+uh(np)^2+sum((uh(1:np-1)-uh(2:np)).^2))/2;

for n = 1 : Ts

    En(n)=0.5*sum(((q0-qm)/k).^2);
    V0(n)=0.5*sum(DV0*q0.*qm);
    qs=Ddiag*q0-qm;
    
    qm=q0;
    q0=qs;
    
    Q(n,:)=q0;
    %w=
    for m = 1 : Nmodes
        w=w+q0(m)*mdShapes(:,:,m);
    end
    W(:,:,n)=w;
    

    
end













%%
figure
plot(t,En-En(3))
hold on
plot(t,V0-V0(3))
figure
E=En+V0;
plot(t,(E-E(3))/E(3))
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
mmoutu=max(max(max(W(:,:,:))));
figure

gif('modeplate.gif') 
for n=2:length(t)-1
surf(x,y,W(:,:,n),"LineWidth",3)
shading interp
zlim([-mmoutu mmoutu])
drawnow
gif
pause(0.001)

end

