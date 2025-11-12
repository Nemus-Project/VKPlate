clear all
close all
clc
%%
%load(['./w10step2dB20modespurpleT=1w0=1e-9f=f0+1.mat'],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')
load(['./w10step2dBpurpleT=1w0=1e-9f=f0+20.mat'],'rho','E','nu','Lz','Lx','Ly','Nmodes','~','Om','~','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','~~','ws','wa')
% load(["./w10step2dB20modespurpleT=2.1w0=1e-8f=f0-30A=?Pa.mat"],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')
% load([ "./w10step2dB20modespurpleT=2.1w0=1e-8f=f0A=100damp=100.mat"],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')
%load(["./chirptest.450modes.mat"],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')

% %%
    % wp=(1/(rho*Lz))*wp;
    % ws=(1/(rho*Lz))*ws;
    % wa=(1/(rho*Lz))*wa;
%load(['./IntrotoQMpurple.mat'],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')


%%

F3=figure('units','normalized','outerposition',[0 0 1 0.7]);

tspec1sp=linspace(0.1,49.5,495)

gif('f0+20.gif')

for iter=1:length(tspec1sp)
Tr=0.5;

tspec1=tspec1sp(iter);
 
tspec2=tspec1+Tr

T1=floor(tspec1*Fs);
T2=floor(tspec2*Fs);

c = linspace(1,10,T2-T1+1);

lim=1.1e-6;

ss = linspace(100,50,T2-T1+1);
turn=linspace(0,360,360);



subplot(1,2,1)
scatter3(wp(T1:T2),ws(T1:T2),wa(T1:T2),"k.");%,s,c,'.')
xlabel("w")
ylabel("dw/dt")
zlabel("d^2w/dt^2")
xlim([-1.5e-6 1.5e-6])
ylim([-2e-3 2e-3])
zlim([-4 4])
set(gca,"FontSize",26)

view(32,42) 

%%% poincarré easy

lspace=Tr*Fs;
testvec=zeros(lspace,1);
%to change the plane, change in 3 places : here change wp to ws
%test if the orbit crosses the wp=0 plane
indtest=1;
for tt=T1:T2
if wp(tt)<0
  if wp(tt+1)>0
      testvec(indtest)=1;
end
end
indtest=indtest+1;
end

% figure
% plot(testvec)

indtest=1;
%compute the points before and after crossing the plane
wpb=zeros(lspace,1);
wsb=zeros(lspace,1);
wab=zeros(lspace,1);
wpa=zeros(lspace,1);
wsa=zeros(lspace,1);
waa=zeros(lspace,1);
for ttt=1:lspace

if  testvec(ttt)==1
    wpb(ttt)=wp(T1+ttt-1);
    wsb(ttt)=ws(T1+ttt-1);
    wab(ttt)=wa(T1+ttt-1);
    wpa(ttt)=wp(T1+ttt);
    wsa(ttt)=ws(T1+ttt);
    waa(ttt)=wa(T1+ttt);
end
end

% finding the intersection with the plane
nbpps=length(nonzeros(waa)); %number of points in the poincarré section

%defining the parametric line that intersect the plane
X0=zeros(nbpps,1);
Y0=zeros(nbpps,1);
Z0=zeros(nbpps,1);

X1=zeros(nbpps,1);
Y1=zeros(nbpps,1);
Z1=zeros(nbpps,1);

itps=1;
for ttt=1:lspace
if  testvec(ttt)==1
    X0(itps)= wpb(ttt);
    Y0(itps)= wsb(ttt);
    Z0(itps)= wab(ttt);

    X1(itps)= wpa(ttt)-wpb(ttt);
    Y1(itps)= wsa(ttt)-wsb(ttt);
    Z1(itps)= waa(ttt)-wab(ttt);
   
  itps=itps+1;  
end
end
%computing the intersection points
testn=[1,0,0]; % change to [0,1,0] another plane of interest
testM=[0,0,0];
ps=zeros(nbpps,3);
for tttt=1:nbpps

testu=[X1(tttt),Y1(tttt),Z1(tttt)];
testN=[X0(tttt),Y0(tttt),Z0(tttt)];

[I,rc]=line_plane_intersection(testu, testN, testn, testM);
ps(tttt,:)=I;
end

%F4=figure%('units','normalized','outerposition',[0 0 1 1]);
testmac=nbpps;
%change the second coordinate of ps to 1
subplot(1,2,2)
scatter(ps(1:testmac,2),ps(1:testmac,3),1000,"k.");%,s,c,'.')
xlabel("dw/dt")
ylabel("d^2w/dt^2")
   xlim([0.1e-3 2e-3] )
  ylim([-2 3.5])
set(gca,"FontSize",26)
drawnow
gif
end