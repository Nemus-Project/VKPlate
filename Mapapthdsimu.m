clear all
close all
clc
%%
% sft1=string(linspace(-10,-1,10));
% sft2=string(["+"+linspace(1,60,60)]);


sft1=string(linspace(-10,-1,10));
sft2=string(["+"+linspace(1,60,60)]);
shft=[sft1,"",sft2];


%shft=["-10","-5","","+5","+10","+15","+20","+25" ];


map=zeros(10,length(shft));
Etot=zeros(10,length(shft)); 
Eover=zeros(10,length(shft)); 
EFund=zeros(10,length(shft));
for iloo=1:length(shft)

    load(strcat("./w10step2dB20modesamberT=0.4w0=1e-9f=f0"+shft(iloo)+"A=0.4.mat"),'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')
   
    if shft(iloo)==""
    f1=Om(1)/(2*pi);
    else
    f1=Om(1)/(2*pi)+str2num(shft(iloo));
    end
    %% time signali
    tplt=linspace(0,T,length(wp));
    wp=(1/(rho*Lz))*wp;
    ws=(1/(rho*Lz))*ws;
    wa=(1/(rho*Lz))*wa;
    %ws=[0;(wp(2:end-1)-wp(1:end-2))*Fs;0];
    %wa=[0;(ws(2:end-1)-ws(1:end-2))*Fs;0];


    %%
    AP=zeros(10,1);
    DIST=zeros(10,1);

    Tr=1;
    f1


    for itit=1:10
        tspec1=2+(itit-1)*5;
        tspec2=tspec1+Tr;
  
        [mod,f]=PSDAM(ws(tspec1*Fs:tspec2*Fs),t(tspec1*Fs:tspec2*Fs),2);

        fs=f(2)-f(1);
        fb=round((f1-10)/fs);
        fBEG=round((f1+10)/fs);
        fEND=round(6000/fs);


        
        Etot(itit,iloo)=sum(mod(1:fEND)); % Fills the Etot matrix with the total psd between 0 and 20kHz
        Eover(itit,iloo)=sum(mod(fBEG:fEND)); % Fills the E matric for each amplitude and frequency
        EFund(itit,iloo)=sum(mod(fb:fBEG));
    end










    for itit=1:10
        tspec1=2+(itit-1)*5;
        tspec2=tspec1+Tr;
        R=yin(ws(tspec1*Fs:tspec2*Fs),Fs);
        ap=R.ap;
        DIST(itit)=thd(ws(tspec1*Fs:tspec2*Fs),Fs);
        AP(itit)=tsnanmean(ap);
    end
    % figure
    % stairs(DIST)
    % title('THD')
    %
    % figure
    % stairs(AP)
    % title('aperiodicity)')

distlm=-13;% -18 for purple, -13.9 for green
APlm=0.005;
    for itit=1:10
        if DIST(itit)> distlm
            map(itit,iloo)=1;
        end
        if AP(itit)>APlm
            map(itit,iloo)=2;
        end
    end
    iloo
end

%% map
load("regimecmap.mat","regimecmap")
nfreq=linspace(-10,60,8);
namp=linspace(1,10,10);
figure('Renderer', 'Painters')
imagesc(nfreq,namp,map)
set(gca,"Ydir","normal","FontSize",26)
xlabel("Normalized frequency (Hz)")
ylabel("Forcing amplitude")
clim([0 2])
colormap(flipud(regimecmap))
colorbar
%%
%load("regimecmap.mat","regimecmap")
nfreq=linspace(-10,60,8);
namp=linspace(1,10,10);
figure('Renderer', 'Painters')
imagesc(nfreq,namp,Etot)
set(gca,"Ydir","normal","FontSize",26)
xlabel("Normalized frequency (Hz)")
ylabel("Forcing amplitude")
 colormap(flipud(bone))
colorbar
%%
%load("regimecmap.mat","regimecmap")
% nfreq=linspace(-10,60,8);
% namp=linspace(1,10,10);
% figure('Renderer', 'Painters')
% imagesc(nfreq,namp,Eover./Etot)
% set(gca,"Ydir","normal","FontSize",26)
% xlabel("Normalized frequency (Hz)")
% ylabel("Forcing amplitude")
% colorbar
%%
% %load("regimecmap.mat","regimecmap")
% nfreq=linspace(-10,60,8);
% namp=linspace(1,10,10);
% figure('Renderer', 'Painters')
% imagesc(nfreq,namp,EFund)
% set(gca,"Ydir","normal","FontSize",26)
% xlabel("Normalized frequency (Hz)")
% ylabel("Forcing amplitude")
% colorbar
%% spectro

load(['./w10step2dB20modesamberT=0.4w0=1e-9f=f0+60A=0.4.mat'],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv','ws','wa')

tplt=linspace(0,T,length(wp));

npwin=2^14;
nlap=npwin/2;
[sv,fmats,tmats,ps] = spectrogram(ws,hann(npwin),nlap,4*npwin,Fs,"ps");




figure
pcolor(tmats,fmats,log10(abs(sv)))

figure('Position', [10 10 1500 500],'Renderer', 'Painters')
pcolor(tmats,fmats(1:10:floor(end/4)),log10(abs(sv(1:10:floor(end/4),:))))

shading interp
colorbar
colormap(flipud(bone))
%caxis([30 100])
%caxis([0.0 0.4])or
%caxis([-3 -1])
ylim([0 4000])
shading interp
colorbar
%caxis([30 100])
%caxis([0.0 0.4])or
%caxis([-3 -1])
%ylim([0 3000])
%clim([0,0.5])
set(gca,"FontSize",26)
ylabel("Frequency (Hz)")
xlabel("Time (s)")
clim([-5 0])


%% check of ap and dist

AP=zeros(10,1);
DIST=zeros(10,1);

Tr=1;

for itit=1:10
    tspec1=2+(itit-1)*5;
    tspec2=tspec1+Tr;
    R=yin(ws(tspec1*Fs:tspec2*Fs),Fs);
    ap=R.ap;
    DIST(itit)=thd(ws(tspec1*Fs:tspec2*Fs),Fs);
    AP(itit)=tsnanmean(ap);
end

figure('Position', [10 10 1500 500],'Renderer', 'Painters')
stairs(DIST,"Linewidth",4,'Color','k')

xlabel("Time (s)")
ylabel("Total Harmonic Distortion")

hold on
yyaxis right
stairs(AP,"Linewidth",4)

ylim([0 0.08])
ylabel("Aperiodicity")

set(gca,"FontSize",26)