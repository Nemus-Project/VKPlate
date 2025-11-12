clear all
close all
clc
addpath '/Users/alexis/Desktop/Thèse Membranes'
%% Loading the data and declaring parameters

Vo=1;

path    = '/Volumes/Backup/LastManips';

%% New detachable setup
SUFF=["100_700Hz60spurplenomousse" ,"100_700Hz60spurplemousse","100_700Hz60sambermoussenomembrane","100_700Hz60sambermousse"];
%SUFF=["100_700Hz60smultimembrane4.3cmsinglepavillon2runsingle2" ,"100_700Hz60smultimembrane4cmsinglepavillon2runsingle2","100_700Hz60smultimembrane3.8cmsinglepavillon2runsingle2","100_700Hz60smultimembrane3.5cmsinglepavillon2runsingle2","100_700Hz60spavillon"];
testhyste=false;%Is this a hysteresis test ?
plotfundout=false;
plotfundin=false;
plotfundv=true;
plotfft=aalse;
plotspectro=false;
plot3sigs=false;


for iter=[length(SUFF),1,2,3]%,3,4]%CHECK THE INDIVIDUAL NUMBERS
suff=SUFF(iter);
tic



filename = strcat(path, '/sweep' ,num2str(Vo), suff, '.mat');%for SweepGreen2
filname=strcat(path, '/sweep' ,num2str(Vo), suff)


load(filename)
Fs=waveGenPar.Fs;
T=waveGenPar.TCapture;
t=time_output;
L=length(t);
f=linspace(100,600,L);


%% Separation of the signals
sigv(:,iter)=sgn_input_avg(:,1)*500;
sigin(:,iter)=sgn_input_avg(:,2)/0.00217;
sigout(:,iter)=sgn_input_avg(:,3)/0.00247;
sigex=sgn_output;
if testhyste==true
    if iter==length(SUFF)
        sigv(:,iter)=flipud(sigv(:,iter));
        sigin(:,iter)=flipud(sigin(:,iter));
        sigout(:,iter)=flipud(sigout(:,iter));
    end
end
%[mod,fpw]=PSDAM(sigout(:,iter),t,30);
[mod,phase,~,fpw]=FFTAM(sigout(:,iter),t,60);

pwx(:,iter)=mod;

npwin=2^15;
nlap=npwin/2;
[s,fmat,tmat,ps] = spectrogram(sigout(:,iter),hann(npwin),nlap,4*npwin,Fs,"ps");

[mod,phase,~,fpw]=FFTAM(sigin(:,iter),t,60);

pwin(:,iter)=mod;

[mod,phase,~,fpw]=FFTAM(sigv(:,iter),t,60);

pwv(:,iter)=mod;

spec(iter,:,:)=s/npwin;



[s,fmat,tmat,ps] = spectrogram(sigin(:,iter),hann(npwin),nlap,4*npwin,Fs,"ps");
specin(iter,:,:)=s/npwin;
%fundin(iter,:)=max(abs(squeeze(specin(iter,:,:))));
[s,fmat,tmat,ps] = spectrogram(sigv(:,iter),hann(npwin),nlap,4*npwin,Fs,"ps");
specv(iter,:,:)=s/npwin;

if iter==length(SUFF)
[fun,ind]=max(abs(squeeze(specin(iter,:,:))));
fundv(iter,:)=fun;
fundin(iter,:)=fun;
fund(iter,:)=max(abs(squeeze(spec(iter,:,:))));
else
    for itee=1:length(ind)
        fund(iter,itee)=abs(spec(iter,ind(itee),itee));
        fundv(iter,itee)=abs(specv(iter,ind(itee),itee));
        fundin(iter,itee)=abs(specin(iter,ind(itee),itee)); 
    end
end

%clear spec
toc
end
%%
% figure
% plot(f,sigout(:,1:end))
% xlabel("Excitation frequency (Hz)")
% ylabel("Microphonic Signal (V)")
% set(gca,"FontSize",26)
% legend("Panneau sans poreux","Poreux épais","Poreux fin")
% %%
% 
% figure
% plot(f,sigin(:,:))
% set(gca,"FontSize",26)
% legend("Pas de panneau","Panneau sans poreux","Poreux épais","Poreux fin")
% 
% %%
% figure
% plot(fpw,pwx(:,1:2))
% hold on
% plot(fpw,pwx(:,3))
% xlabel("Frequency (Hz)")
% ylabel("Power")
% xlim([100 2500])
% set(gca,"FontSize",26)
% legend("Panneau seul","Green seul","Poreux+Purple")
% 
% 
% %%
% 
% %%
% figure
% plot(fpw,pwin(:,1:end))
% xlabel("Frequency (Hz)")
% ylabel("Power")
% xlim([100 2500])
% set(gca,"FontSize",26)
% legend("Panneau seul","Purple seul","Poreux+Purple")



%% panneau poreux
% newcolors = [
%              0.8 0.4 0
%              0.69 0.28 0.82
%              0 0.7 0
%              0 0.4470 0.7410
%              0 0 0
%              0.1 0.50 0.95
%              0 0 0 ];%multisize

newcolors = [
             0.75 0.75 0.75
             0.85   0.64  0.12
             0.8 0.4 0
             0.69 0.28 0.82
             0 0.4470 0.7410
             0 0 0
             0.1 0.50 0.95
             0 0 0 ];%newmembranes
if testhyste==true
    newcolors = [
             1 0 1
             0.5 0 0.5
             0 0 0
             0.1 0.50 0.95
             0 0 0 ]; % Purple hyste

     % newcolors = [
     %         0 1 0
     %         0 0.5 0
     %         0 0 0
     %         0.1 0.50 0.95
     %         0 0 0 ]; %Green hyste
      % 
      % newcolors = [
      %        1 0.64 0
      %        0.8 0.5 0
      %        0 0 0
      %        0.1 0.50 0.95
      %        0 0 0 ]; %Green hyste
end
%% FUND
if plotfundout==true
figure
colororder(newcolors)
ffund = linspace(100,700,length(fund));
plot(ffund,20*log10(fund(1:end-1,:)/2e-5),"LineWidth",4)
hold on
plot(ffund,20*log10(fund(end,:)/2e-5),"LineWidth",4,"Color","k","LineStyle","-.")
xlabel("Frequency (Hz)")
ylabel("Power(dB)")
xlim([100 700])
set(gca,"FontSize",26)
legend("Tan","Purple","Green","Blue","NoMembrane")
if testhyste==true
    legend("Purple forward","Purple backward","No Membrane")
end
end
%% FUNDIN
if plotfundin==true
figure
colororder(newcolors)
ffund = linspace(100,700,length(fund));
plot(ffund,20*log10(fundin(1:end-1,:)/2e-5),"LineWidth",4)
hold on
plot(ffund,20*log10(fundin(end,:)/2e-5),"LineWidth",4,"Color","k","LineStyle","-.")
xlabel("Frequency (Hz)")
ylabel("Power(dB)")
xlim([100 700])
set(gca,"FontSize",26)
legend("Tan","Purple","Green","Blue","NoMembrane")
if testhyste==true
    legend("Purple forward","Purple backward","No Membrane")
end
end
%% FUNDV multiple thicknesses
% 
% if plotfundv==true
% figure
% colororder(newcolors)
% ffund = linspace(100,700,length(fund));
% plot(ffund,(fundv(1:end-1,:)),"LineWidth",4)
% xlabel("Frequency (Hz)")
% ylabel("Amplitude (mm.s^{-1})")
% xlim([100 700])
% set(gca,"FontSize",26)
% legend("Tan","Purple","Green","Blue","NoMembrane")
% %legend("4.3cm","4cm","3.8cm","3.5cm","NoMembrane")
% title("Fundamental of the vibration")
% %xlim([100 500])
% if testhyste==true
%     legend("Purple forward","Purple backward","No Membrane")
% end
% end

%% FUNDV multiple

newcol2=[0.8 0.4 0
         0 0 0
         0.8 0.4 0
         0 0 0
         0 0 0
         0.1 0.50 0.95
         0 0 0 ];
if plotfundv==true
figure
colororder(newcolors)
ffund = linspace(100,700,length(fund));
plot(ffund,(fundv(1,:)),"LineWidth",4)
hold on
plot(ffund,(fundv(2,:)),"LineWidth",4)
plot(ffund,(fundv(3,:)),"LineWidth",4)%,"LineStyle","--")
plot(ffund,(fundv(4,:)),"LineWidth",4)%,"LineStyle","--")
xlabel("Frequency (Hz)")
ylabel("Amplitude (mm.s^{-1})")
xlim([100 700])
set(gca,"FontSize",26)
%legend("Tan","Purple","Green","Blue","NoMembrane")
legend("10\mum","20\mum","25\mum","40\mum","NoMembrane")
title("Fundamental of the vibration")
xlim([100 500])
if testhyste==true
    legend("Purple forward","Purple backward","No Membrane")
end
end

%% FFT
if plotfft==true
figure
colororder(newcolors)
plot(fpw,pwv(:,1:end),"LineWidth",2)
xlabel("Frequency (Hz)")
ylabel("Power")
xlim([100 5000])
set(gca,"FontSize",26)
legend("Purple","Purple+Poreux","NoMembrane")
if testhyste==true
    legend("Purple forward","Purple backward","No Membrane")
end
end
%% SPECTRO
if plotspectro==true
% spectro
ite=1;
for nst=[2];
steeeeeee=abs(squeeze(specv(nst,:,:)));
%stest=20*log10(2*steeeeeee/(2e-5));
stest=steeeeeee;
figure
subplot(3,1,[1,2])
pcolor(tmat,fmat,stest)

shading interp
colorbar
%caxis([30 100])
%caxis([0.0 0.4])or
%caxis([-3 -1])
ylim([0 5000])
set(gca,"FontSize",26)
ylabel("Frequency (Hz)")
xlabel("Time (s)")
%set (gca,'Xdir','reverse')
title ("Sweep from 100 to 700 Hz micin multimembrane")
colormap(flipud(bone))

subplot(3,1,3)
plot(t,sigout(:,nst))
set(gca,"FontSize",26)
ite=ite+1;
end
end
%% Signals
if plot3sigs==true
ite=1;
figure
for nst=1:length(SUFF);
subplot(length(SUFF),1,ite)
%plot(t,abs(sigout(:,nst))+abs(sigin(:,nst)),t,sigin(:,nst),t,sigout(:,nst))
plot(t,sigin(:,nst),t,sigout(:,nst))
hold on
ylim([-10 10])
set(gca,"FontSize",26)
xlabel("Time (s)")
%ylabel("Pressure (Pa)")
%title("Mic Out")
ite=ite+1;
end
end
% %%
% for nst=1:length(SUFF)
% integ(nst)=trapz(abs(sigout(:,nst)))
% integin(nst)=trapz(abs(sigin(:,nst)))
% 
% end

%% poreux
figure
%colororder(newcolors)
ffund = linspace(100,700,length(fund));
plot(ffund,(fund(1,:)),"LineWidth",4)
hold on
plot(ffund,(fund(2,:)),"LineWidth",4)
plot(ffund,(fund(3,:)),"LineWidth",4)
legend("10\mum","10\mum foam","Only foam","40\mum","NoMembrane")
set(gca,"FontSize",26)