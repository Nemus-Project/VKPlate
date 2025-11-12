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
%% time signal
tplt=linspace(0,T,length(wp));

figure
plot(tplt,wp)
% ws=[0;(wp(2:end-1)-wp(1:end-2))*Fs;0];
% wa=[0;(ws(2:end-1)-ws(1:end-2))*Fs;0];
%fplt=linspace(100,2500,44100*30);
%%
figure
plot(tplt(1:10:end),ws(1:10:end))
%xlim([100 2500])
xlabel("Time (s)")
ylabel("Velocity (m.s^-1)")
set(gca,'FontSize',24)


%% envelopes of modes
% [yupper,ylower] = envelope(Qsv(:,:),8000,'analytic');
%  figure
% plot(tplt,yupper)
%% Spectrogram velocity

% [y,fs]=audioread("Intro_to_QM-2.mp3");
% yy=y(1:fs*19,1);

npwin=2^12;
nlap=npwin/2;
[sv,fmats,tmats,ps] = spectrogram(ws,hann(npwin),nlap,4*npwin,Fs,"ps");




figure%('Position', [10 10 1500 500],'Renderer', 'Painters')
pcolor(tmats,fmats(1:5:floor(end/2)),log10(abs(sv(1:5:floor(end/2),:))))

shading interp
colorbar
colormap(flipud(bone))
%caxis([30 100])
%caxis([0.0 0.4])or
%caxis([-3 -1])
ylim([0 4000])
%clim([0,0.5])
set(gca,"FontSize",26)
ylabel("Frequency (Hz)")
xlabel("Time (s)")
%set (gca,'Xdir','reverse')

%colormap(flipud(bone))

clim([-4 1])
[fun,ind]=max(abs(squeeze(s(:,:))));
ind2=linspace(ind(1),ind(end),length(ind));
fund=zeros(length(ind),1);
%%

tind=5;
tnum= ceil(tind*length(tmats)/50);
figure
%plot(fmat(1:10:floor(end/4)),log10(abs(s(1:10:floor(end/4),tnum))),"Color",[0.69 0.28 0.82],"LineWidth",2)
plot(fmats,log10(abs(sv(:,tnum))),"Color","b","LineWidth",2)
xlim([10 5000])
%ylim([0 2.5])
set(gca,"FontSize",26)
xlabel("Frequency (Hz)")
ylabel("Amplitude")


%%
[fun,ind]=max(abs(squeeze(s(:,:))));
ind2=linspace(ind(1),ind(end),length(ind));
fund=zeros(length(ind),1);

%%
Tr=1;

tspec1=0.1;
 
tspec2=tspec1+Tr
c = linspace(1,10,tspec2*Fs-tspec1*Fs+1);

lim=1.1e-6;

ss = linspace(100,50,tspec2*Fs-tspec1*Fs+1);
turn=linspace(0,360,360);


F3=figure%('units','normalized','outerposition',[0 0 1 1]);

scatter3(wp(tspec1*Fs:tspec2*Fs),ws(tspec1*Fs:tspec2*Fs),wa(tspec1*Fs:tspec2*Fs),"k.");%,s,c,'.')
xlabel("w")
ylabel("dw/dt")
zlabel("d^2w/dt^2")
% xlim([-1.2e-6 1.2e-6])
% ylim([-2e-3 2e-3]);
% zlim([-5 5]);
% hold on
% patch([0 0 0 0], [1.5e-3 1.5e-3 0 0 ], [3 -3 -3 3], [0.8,0.8,0.8],'FaceAlpha',.85)  



%view(70.5,5.3947)
% gif('./img/quasip.gif');
% for tu=1:length(turn)
%view(turn(tu),5.3947)
% gif
set(gca,"FontSize",26)
% end
%view(70,5.3947) % for periodic

view(24,67) % for aperiodic
% hold on
% 
% plot3(T0OrP(1,:),T0OrP(2,:),T0OrP(3,:))
% plot3(U0OrP(1,:),U0OrP(2,:),U0OrP(3,:))
% plot3(V0OrP(1,:),V0OrP(2,:),V0OrP(3,:))
% legend("T0","U0","V0")
%colormap(Transi)


% poincarré section :
% xp=linspace(0,1);
% yp=linspace(0,1);
% zp=linspace(0,1);
% ps=tspec1*Fs; %point of the section
% T0=[prooot(ps,1),prooot(ps,2),prooot(ps,3)]-[prooot(ps+1,1),prooot(ps+1,2),prooot(ps+1,3)] %tangent vector
% T0P=T0.*linspace(0,1)';
% U0=[1,0,0]*1e-6;
% U0P=U0.*linspace(0,1)';
% V0=[0,1,0]*1e-6;
% V0P=V0.*linspace(0,1)';
% figure
% plot3(T0P(:,1),T0P(:,2),T0P(:,3))
% hold on
% plot3(U0P(:,1),U0P(:,2),U0P(:,3))
% plot3(V0P(:,1),V0P(:,2),V0P(:,3))
% legend("T0","U0","V0")
% 
% BAS=[T0',U0',V0'];
% [BASorth,OrtN]=GramSchmidt(BAS);
% T0Or=BASorth(:,1);
% U0Or=BASorth(:,2);
% V0Or=BASorth(:,3);
% 
% p1=prooot(ps,1);
% p2=prooot(ps,2);
% p3=prooot(ps,3);
% 
% l1=linspace(p1,p1*1e-6);
% l2=linspace(p2,p2*1e-6);
% l3=linspace(p3,p3*1e-6);
% %%

% T0OrP=[T0Or(1).*l1;T0Or(2).*l2;T0Or(3).*l3];
% U0OrP=[U0Or(1).*l1;U0Or(2).*l2;U0Or(3).*l3];
% V0OrP=[V0Or(1).*l1;V0Or(2).*l2;V0Or(3).*l3];
% 
% figure
% plot3(T0OrP(1,:),T0OrP(2,:),T0OrP(3,:),"LineWidth",4)
% hold on
% plot3(U0OrP(1,:),U0OrP(2,:),U0OrP(3,:),"LineWidth",4)
% plot3(V0OrP(1,:),V0OrP(2,:),V0OrP(3,:),"LineWidth",4)
% legend("T0","U0","V0")
% 
% %% Intersection to the plane U0 V0
% [x y] = meshgrid(-1e-6:1e-7:1e-6);
% z = (T0(1).*x + T0(2).*y)/T0(3);
% figure
% surf(x,y,z);

%%% poincarré easy


%Tr=3;

%tspec1=21;
 
%tspec2=tspec1+Tr;

lspace=Tr*Fs;
testvec=zeros(lspace,1);
%to change the plane, change in 3 places : here change wp to ws
%test if the orbit crosses the wp=0 plane
indtest=1;
for tt=tspec1*Fs:tspec2*Fs
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
    wpb(ttt)=wp(tspec1*Fs+ttt-1);
    wsb(ttt)=ws(tspec1*Fs+ttt-1);
    wab(ttt)=wa(tspec1*Fs+ttt-1);
    wpa(ttt)=wp(tspec1*Fs+ttt);
    wsa(ttt)=ws(tspec1*Fs+ttt);
    waa(ttt)=wa(tspec1*Fs+ttt);
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

F4=figure%('units','normalized','outerposition',[0 0 1 1]);
testmac=nbpps;
%change the second coordinate of ps to 1
scatter(ps(1:testmac,2),ps(1:testmac,3),1000,"k.");%,s,c,'.')
xlabel("dw/dt")
ylabel("d^2w/dt^2")
%   xlim([0.5e-3 1.5e-3] )
 % ylim([-2 3])
set(gca,"FontSize",26)


%% manifold

Tr=2;

redprot=zeros(6*Tr*Fs,3);
F3=figure('Renderer','Painter') %('units','normalized','outerposition',[0 0 1 1]);
hold on
for itit=[1,2,3,4,8,9,10]
tspec1=2+(itit-1)*5;
 
tspec2=tspec1+Tr;

redprot(1+(itit-1)*Tr*Fs:(itit)*Tr*Fs+1,:)=[wp(tspec1*Fs:tspec2*Fs),ws(tspec1*Fs:tspec2*Fs),wa(tspec1*Fs:tspec2*Fs)];
if itit==3
    scatter3(wp(tspec1*Fs:tspec2*Fs),ws(tspec1*Fs:tspec2*Fs),wa(tspec1*Fs:tspec2*Fs),300,'k.')
elseif itit==5
    scatter3(wp(tspec1*Fs:tspec2*Fs),ws(tspec1*Fs:tspec2*Fs),wa(tspec1*Fs:tspec2*Fs),300,'k.')

elseif itit==9
    scatter3(wp(tspec1*Fs:tspec2*Fs),ws(tspec1*Fs:tspec2*Fs),wa(tspec1*Fs:tspec2*Fs),300,'k.')

else
scatter3(wp(tspec1*Fs:tspec2*Fs),ws(tspec1*Fs:tspec2*Fs),wa(tspec1*Fs:tspec2*Fs),300,'.')
end
xlabel("w(t)")
ylabel("dw/dt")
zlabel("dw^2/dt^2")
view(70.5,5.3947)
set(gca,"FontSize",26)
end

%%
% figure
% plot3(redprot(:,1),redprot(:,2),redprot(:,3),'k.')
% 
% rprotx=redprot(:,1);
% rproty=redprot(:,2);
% rprotz=redprot(:,3);
% 
% [rprotxq,rprotyq] = meshgrid(min(rprotx):1e-9:max(rprotx),min(rproty):1e-5:max(rproty));
% rprotzq = griddata(rprotx,rproty,rprotz,rprotxq,rprotyq);
% figure
% mesh(rprotxq,rprotyq,rprotzq);
% %colormap(UniPur)
% hold on
% plot3(redprot(:,1),redprot(:,2),redprot(:,3),'k')



%%
%  figure
%  plot3(redprot(:,1),redprot(:,2),redprot(:,3),'k.')
% view(70.5,5.3947)
%%
fspec=128.5133+15
etotv=zeros(1,10);
    eoverv=zeros(1,10);
    eFundv=zeros(1,10);
for n =0:9 % This loops through the ten steps of the excitation and windows each step
   
    lowt=Fs*(n*5+1); %begining of the time window one sec after the change of excitation
    hight=Fs*(n*5+3); % end of the time window two secs after the begining
    
    
    [modv,f]=PSDAM(ws(lowt:hight),t(lowt:hight),2);   %power spectral density of the vibration
    

    
    fs=f(2)-f(1);
    nHarm=round(20000/fspec);
    fb=round((fspec-20)/fs);
    fBEG=round((fspec+20)/fs);
    %fBEG=round((500)/fs);
    fharm=linspace(fspec*2,fspec*nHarm,nHarm-1);
    fEND=round(20000/fs); %end frequency is 20 kHz
   %R=yin(sigout(lowt:hight),P);
    
    
  

    etotv(n+1)=sum(modv(1:fEND));
    eoverv(n+1)=sum(modv(fBEG:fEND));
    eFundv(n+1)=sum(modv(fb:fBEG));    
end


figure
plot(etotv,LineWidth=4)
hold on
plot(eFundv,LineWidth=4)
plot(eoverv,LineWidth=4)
legend('total','fund','over')
function [Q, R] = GramSchmidt(X)
% Modified Gram-Schmidt orthonormalization (numerical stable version of Gram-Schmidt algorithm) 
% which produces the same result as [Q,R]=qr(X,0)
% Written by Mo Chen (sth4nth@gmail.com).
[d,n] = size(X);
m = min(d,n);
R = zeros(m,n);
Q = zeros(d,m);
for i = 1:m
    v = X(:,i);
    for j = 1:i-1
        R(j,i) = Q(:,j)'*v;
        v = v-R(j,i)*Q(:,j);
    end
    R(i,i) = norm(v);
    Q(:,i) = v/R(i,i);
end
R(:,m+1:n) = Q'*X(:,m+1:n);
end
