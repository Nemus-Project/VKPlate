clear all 
close all
clc
% This code solves the Duffing Oscillator with the quadratization of enery method and plots
% , we can see that the energy is conserved, while in the first code it was not for some mysterious reason 

%THIS CODE WORKS, DONT MODIFY IT
%% Declaration of parameters
omeg0=10; %Frequency of the oscilator
T=10; %Final time in seconds
k=0.00001; %Timestep
t=linspace(0,T,floor(T/k)+1); %Time vector
th=linspace(k/2,T+k/2,floor(T/k)+1); %Interleaved time vector
Hsav=zeros(1,length(t)); %Leapfrog discrete energy
vlin=zeros(1,length(t)); %Leapfrog discrete energy
u=zeros(1,length(t)); %Displacement vector
psi=zeros(1,length(t)); %Interleaved displacement vector
g=zeros(1,length(t));
gamma=0;
A=1;
um1=A;
[s,c,d] = ellipj(sqrt(omeg0^2+gamma*A^2)*t,(gamma*A^2)/(2*gamma*A^2+2*omeg0^2));
uana=A*c;
u0=uana(2);
[s,c12,d] = ellipj(sqrt(omeg0^2+gamma*A^2)*(k/2),(gamma*A^2)/(2*gamma*A^2+2*omeg0^2));
u12=A*c12;
V12=(gamma*u12^4)/4;
h0=((u0-um1)^2)/(2*k^2)+(omeg0^2*u12^2)/2+V12; 
eta=1.919e-2;
eps=eta*h0;
u(1)=um1;
u(2)=u0;
up1=0;
psim1=sqrt(2*V12+eps);
psi(1)=psim1;


%% Fast method
for n=2:length(t)-1
    %u(n+1)=u(n)*(2-k^2*omeg0^2)-u(n-1);
    Hsav(n-1)=(1/2)*((u0-um1)/k)^2+(omeg0^2/2)*um1*u0+(1/2)*(psim1)^2;
    %vlin(n-1)=;
    g0=(gamma*u0^3)/(sqrt(gamma*(u0^4)/2+eps));
    %up1=(1/(1+(g0^2*k^2)/4))*(2*u0-um1+((g0^2*k^2)/4)*um1-g0*k^2*psim1);%add the linear term, maybe even the exact integration
    %up1=2*cos(omeg0*k)*u0-um1-k^2*g0*psim1;
    %up1=(2-k^2*omeg0^2)*u0-um1-k^2*g0*psim1;
    %up1=(1/(1+(g0^2*k^2)/4))*((2-k^2*omeg0^2)*u0+(((g0^2*k^2)/4)-1)*um1-k^2*g0*psim1);
    up1=(1/(1+((omeg0^2*k^2)/2)+(g0^2*k^2)/4))*(2*u0-(1+((omeg0^2*k^2)/2)-(g0^2*k^2)/4)*um1-k^2*g0*psim1);
    
    psi0=psim1+(g0/2)*(-um1+up1);
    
    g(n)=g0;
    psim1=psi0;
    psi(n)=psi0;
    um1=u0;
    u0=up1;
    u(n+1)=up1;
end

%%
figure
plot(t,u,LineWidth=3)
hold on
plot(t,uana,LineWidth=3,LineStyle="--")
%plot(t,g,LineWidth=3)
legend("u","uana")
xlabel("Time(s)")
ylabel("Amplitude")
set(gca,'Fontsize',20)
%xlim([0,10])
%ylim([-10,10])
%%
% figure
% plot(u(3:length(t)-1),u(4:length(t)))

%%
%%
figure
plot(th(1:1000),(Hsav(1:1000)-Hsav(3))/Hsav(3),LineWidth=3)
