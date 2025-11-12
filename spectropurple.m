clear all
%close all
clc


%%
load('/Volumes/BolognaDisk/2024/Steps/05V/sweep2100_2500Hz30spurplerun10.mat')
Fs=waveGenPar.Fs;
T=waveGenPar.TCapture;
t=time_output;
L=length(t);



sigv=sgn_input_avg(:,1);
sigin=sgn_input_avg(:,2);
sigout=sgn_input_avg(:,3);
sigex=sgn_output;
sigout=sigout;
npwin=2^13;
nlap=npwin/4;
[s,fmat,tmat] = spectrogram(sigv,hann(npwin),nlap,4*npwin,Fs); %CHANGE SIGNAL HERE
nm=length(s(1,:));
nn=length(s(:,1));
ymat=linspace(0,T,nm);
[Xmat,Ymat]=meshgrid(ymat,fmat);
%% Plots

freqplt=linspace(100,2500,Fs*30);
figure
plot(freqplt,sigv)
xlim([100 2500])
set(gca,"FontSize",26)
%%
figure%('Renderer','Painter')
%subplot(3,2,1:4)
pcolor(tmat,fmat,log10(abs(s)))
%surf(tmat,fmat,log10(abs(s)))
shading interp
colorbar
caxis([-2 1])
ylim([0 4000])
set(gca,"FontSize",26)
ylabel("Frequency (Hz)")
xlabel("Time (s)")
%title (strcat("Spectrogram of the ", num2str(f1)," Hz signal"))
colormap(flipud(bone))