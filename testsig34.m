clear all
%close all
clc
%% Loading the data and declaring parameters
%load('PdB05V.mat')

Vo=0.05;

% mem='multimembranenomembrane';
% suff=['Hz50s0.5Vrunpavillon3'];

f1=150;
    toc
    tic
%f1=230;

% path    = '/Volumes/BolognaDisk/2024/Steps/05V';


path    = '/Volumes/Backup/LastManips';

%mem='ambermoussenomembrane';
mem='multimemnomoussenomembrane';
suff='Hz50s0.5Vrunjuly';


filename = [path '/10step' [mem  num2str(f1) suff] '.mat'];%for SweepGreen2
filname=[path '/step0.05V' num2str(f1) 'Hz50s']


load(filename)
Fs=waveGenPar.Fs;
T=waveGenPar.TCapture;
t=time_output;
L=length(t);
Fmin=500;
Fmax=f1;
%Fmax=f1;
sigoutnomem=sgn_input_avg(:,3);
%sig=sigout(lowt:hight)*10/0.091;
%% sig memnomousse

mem='multimempurplenomousse';
filename = [path '/10step' [mem  num2str(f1) suff] '.mat'];%for SweepGreen2
load(filename)

sigv=sgn_input_avg(:,1);
sigin=sgn_input_avg(:,2);
sigout=sgn_input_avg(:,3);
sigex=sgn_output;
sigout=sigout;
npwin=2^15;
nlap=npwin/4;
[s,fmat,tmat] = spectrogram(sigout,hann(npwin),nlap,4*npwin,Fs); %CHANGE SIGNAL HERE
nm=length(s(1,:));
nn=length(s(:,1));
ymat=linspace(0,T,nm);
[Xmat,Ymat]=meshgrid(ymat,fmat);
%% sig moussesolo
mem='multimemmoussenomembrane';
filename = [path '/10step' [mem  num2str(f1) suff] '.mat'];%for SweepGreen2
load(filename)


sigoutms=sgn_input_avg(:,3);


%% sig mem+mousse

mem='multimempurplemousse';
filename = [path '/10step' [mem  num2str(f1) suff] '.mat'];%for SweepGreen2
load(filename)


sigoutmm=sgn_input_avg(:,3);



%%
% Thdis=thd(sigout,Fs,15);
% 
%
env=envelope(sigex);
figure
plot(t(5000:50:end-7000),env(5000:50:end-7000),"LineWidth",8,"Color","k")
set(gca,"FontSize",26)


%% YIN tests


SIG=zeros(Fs*2+1,10);
for n=0:9
    tic
lowt=Fs*(n*5+1);
hight=Fs*(n*5+3);
sig=sigout(lowt:hight)*10/0.091;
sigdB=20*log10(abs(sig)/2e-5);
s_sqr=(abs(sigout)).^2;
Es=sum(s_sqr);
PdB(n+1)=20*log10(rms(sig)/2e-5);
SIG(:,n+1)=sigout(lowt:hight);
SIGNM(:,n+1)=sigoutnomem(lowt:hight);
SIGMS(:,n+1)=sigoutms(lowt:hight);
SIGMM(:,n+1)=sigoutmm(lowt:hight);
SIGEX(:,n+1)=sigex(lowt:hight);
 SIGFILT(:,n+1)=bandpass(sigout(lowt:hight),[f1-5, f1+5],Fs);
 SIGFILT2(:,n+1)=bandpass(sigout(lowt:hight),[f1+50, 20000],Fs);
 SIGH2(:,n+1)=bandpass(sigout(lowt:hight),[2*f1-5, 2*f1+5],Fs);
 SIGH3(:,n+1)=bandpass(sigout(lowt:hight),[3*f1-5, 3*f1+5],Fs);
 SIGH4(:,n+1)=bandpass(sigout(lowt:hight),[4*f1-5, 4*f1+5],Fs);

Thdis(:,n+1)=thd(sigout(lowt:hight),Fs,15)
toc
% SIGFILT3(:,n+1)=bandpass(sigout(lowt:hight),[f1-5, 20000],Fs);
end
PdB
%% THD
% figure
% stairs([Thdis,Thdis(end)],"LineWidth",3)


%%
 % Mex=max(SIGEX(:,:));
 % Mnm=max(SIGNM(:,:));
 % Msig=max(SIG(:,:));

Mex=sqrt(mean(SIGEX(:,:).^2));
Mnm=sqrt(mean(SIGNM(:,:).^2));
Msig=sqrt(mean(SIG(:,:).^2));
MsigMS=sqrt(mean(SIGMS(:,:).^2));
MsigMM=sqrt(mean(SIGMM(:,:).^2));

 Mfilt=sqrt(mean(SIGFILT(:,:).^2));
 Mfilt2=sqrt(mean(SIGFILT2(:,:).^2));
 MH2=sqrt(mean(SIGH2(:,:).^2));
MH3=sqrt(mean(SIGH3(:,:).^2));
MH4=sqrt(mean(SIGH4(:,:).^2));

% Mfilt3=sqrt(mean(SIGFILT2(:,:).^2))

%% THD and harmonics plot
% 
% tplt=linspace(0,50,10);
% figure('Position', [10 10 2000 500])
% plot(tplt,Mfilt,"LineWidth",3,"Marker","o","MarkerSize",15)
% hold on
% plot(tplt,MH2,"LineWidth",3,"Marker","o","MarkerSize",15)
% plot(tplt,MH3,"LineWidth",3,"Marker","o","MarkerSize",15)
% %plot(tplt,MH4,"LineWidth",3,"Marker","o","MarkerSize",15)
% yyaxis right
% plot(tplt,Thdis,"LineWidth",3,"Color","k","Marker","diamond","MarkerSize",15)
% 
% legend("Fundamental","Harmonic 2",'Harmonic 3',"THD")
% 
% set(gca,"FontSize",26)
%% Multiple Signals
nint=4; % specific amplitude of interest
tp=t(1:2*Fs+1)*f1;
Fivep=find(floor(tp)==5,1);%index of the 5th period

ofNM=mean(SIGNM(1:Fivep,nint)); %Offsets
of=mean(SIG(1:Fivep,nint));
ofMS=mean(SIGMS(1:Fivep,nint));
ofMM=mean(SIGMM(1:Fivep,nint));

SIGNM=SIGNM-ofNM;
SIG=SIG-of;
SIGMS=SIGMS-ofMS;
SIGMM=SIGMM-ofMM;

% Signaux individuels

figure
 tprig=tp;
p=plot(tp,SIGNM(:,[nint])/10,'--',tprig,SIG(:,nint)/10,'-');
hold on
plot(tprig,SIGMS(:,nint)/10,'-',"LineWidth",3,"Color",[0.349, 0.522, 0.392])
plot(tprig,SIGMM(:,nint)/10,'-',"LineWidth",3,"Color","r")
p(1).LineWidth = 3;
p(1).Color = [0    0    0];
p(2).LineWidth = 2;
p(2).Color = [0.69 0.28 0.82];

set(gca,"FontSize",26)
ylabel("Microphonic signal (V)")
xlabel("Number of periods")
xlim([0 5])
legend('Reference','Membranes','Foam',"Membranes+Foam")
%ylim([-3e-3 3e-3])
%ylim([-12e-3 12e-3])
ylim([-10e-3 10e-3])
%title (strcat("Signal at ", num2str(f1), "Hz, ", mem, " plate" ))


%%
% figure
% p=plot(tp,SIGEX(:,[1,6,10]));
% 
% p(1).LineWidth = 2;
% p(2).LineWidth = 3;
% p(3).LineWidth = 4;
% ylabel("Excitation signal (V)")
% xlabel("Number of periods")
% xlim([0 5])
%%

% figure
% p=plot(tp,SIGFILT(:,[1,6,10]));
% 
% p(1).LineWidth = 2;
% p(2).LineWidth = 3;
% p(3).LineWidth = 4;
% ylabel("Excitation signal (V)")
% xlabel("Number of periods")
% xlim([200 205])
%% FITS
nfits=10;
 
 fitnm = polyfit(Mex,Mnm,1)
 fitnmplt =polyval(fitnm,Mex)
 
x1sig=linspace(Mex(1),Mex(10),3)
fitsig= polyfit(Mex(1:nfits),Msig(1:nfits),1)
fitsigplt =polyval(fitsig,x1sig)

 figure('Position', [10 10 2000 500])

       
 %scatter(Mex,Mnm,150,"k","Marker","o","LineWidth",3)

 semilogx(Mex,20*log10(Mnm*(10/0.091)/(2e-5)/10),"k","Marker","o","LineWidth",4,"MarkerSize",20)
 hold on
 %plot(Mex,fitnmplt,"Color","k","LineStyle","--","LineWidth",3)
 %scatter(Mex,Msig,150,"r","Marker","diamond","LineWidth",3)
semilogx(Mex,20*log10(Msig*(10/0.091)/(2e-5)/10),"Color",[0.69 0.28 0.82],"Marker","^","LineWidth",4,"MarkerSize",20)
semilogx(Mex,20*log10(MsigMS*(10/0.091)/(2e-5)/10),"Color",[0.349, 0.522, 0.392],"Marker","v","LineWidth",4,"MarkerSize",20)
semilogx(Mex,20*log10(MsigMM*(10/0.091)/(2e-5)/10),"r","Marker","diamond","LineWidth",4,"MarkerSize",20)

%semilogx(Mex,20*log10(Mfilt*(10/0.091)/(2e-5)),"Color","#E96703","Marker","x","LineWidth",4,"LineStyle","--","MarkerSize",20)
%semilogx(Mex,20*log10(Mfilt2*(10/0.091)/(2e-5)),"Color","#C4526C","Marker","x","LineWidth",4,"LineStyle","--","MarkerSize",20)
%plot(x1sig,fitsigplt,"Color","r","LineStyle","--","LineWidth",3)
%scatter(Mex,Mfilt3,150,"k","Marker","o","LineWidth",2)
%scatter(Mex,Mfilt,150,"b","Marker","diamond","LineWidth",2)
ylabel("RMS of the microphonic signal (V)")
xlabel("RMS of the input signal (V)")
legend('Reference','Membranes','Foam',"Membranes+Foam")
title (strcat("RMS of the signals at ", num2str(f1), "Hz, ", mem, " plate" ))
%legend("RMS without plate",['RMS ' mem ' plate'],"RMS at the fundamental","RMS over the fundamental")
set(gca,"FontSize",26)
%%

figure('Position', [10 10 2000 500])

 semilogx(Mex,20*log10(Mfilt*(10/0.091)/(2e-5)/10),"b","Marker","square","LineWidth",4,"MarkerSize",20)
 hold on

 semilogx(Mex,20*log10(Mfilt2*(10/0.091)/(2e-5)/10),"g","Marker","+","LineWidth",4,"MarkerSize",20)

ylabel("RMS of the microphonic signal (V)")
xlabel("RMS of the input signal (V)")

title (strcat("RMS of the signals at ", num2str(f1), "Hz, ", mem, " plate" ))
legend("Fundamental RMS","Higher than fund RMS")
set(gca,"FontSize",26)
%%
