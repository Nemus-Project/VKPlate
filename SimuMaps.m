clear all
close all
clc
%%
load('./w10step2dBurpleT=1w0=1e-9f=f0.mat','rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv','Npsi','wp','Fs','T','t','s','fmat','tmat','faxis','mod','A','Qsv')


%% time signal
tplt=linspace(0,T,length(wp));

figure
plot(tplt,wp)
%% envelopes of modes
% [yupper,ylower] = envelope(Qsv(:,:),8000,'analytic');
% figure
% plot(tplt,yupper)


%% Spectrogram



figure
pcolor(tmat,fmat,log10(abs(s)))

shading interp
colorbar
clim([-10 -3.5])
ylim([0 4000])

set(gca,"FontSize",26)
ylabel("Frequency (Hz)")
xlabel("Time (s)")
%set (gca,'Xdir','reverse')
title ("Spectro SV")

[fun,ind]=max(abs(squeeze(s(:,:))));
ind2=linspace(ind(1),ind(end),length(ind));
fund=zeros(length(ind),1);
%fund(iter,:)=max(abs(squeeze(spec(iter,:,:))));
% 
% for itee=1:length(ind)
%         fund(itee)=abs(s(floor(ind2(itee)),itee));
% end
% figure
% plot(fund,"LineWidth",3)
% xlabel("Time")
% ylabel("Amplitude")