clear all 
close all
clc
% This code solves the Duffing Oscillator with the quadratization of enery method and plots
% , we can see that the energy is conserved, while in the first code it was not for some mysterious reason 

%THIS CODE WORKS, DONT MODIFY IT
%% Declaration of parameters
omeg0=10;                                                                           %Frequency of the oscilator
T=10;                                                                               %Final time in seconds
k=0.001;                                                                            %Timestep
t=linspace(0,T,floor(T/k)+1);                                                       %Time vector
th=linspace(k/2,T+k/2,floor(T/k)+1);                                                %Interleaved time vector
Hsav=zeros(1,length(t));                                                            %Discrete energy vector
vlin=zeros(1,length(t));                                                            %Linear potential energy
klin=zeros(1,length(t));                                                            %Kinetic energy vector
u=zeros(1,length(t));                                                               %Displacement vector
psi=zeros(1,length(t));                                                             %Psi vector
g=zeros(1,length(t));                                                               %G vector
gamma=300;                                                                          %NL Parameter
A=1;                                                                                %Amplitude of initial displacement
um1=A;
[s,c,d] = ellipj(sqrt(omeg0^2+gamma*A^2)*t,(gamma*A^2)/(2*gamma*A^2+2*omeg0^2));    %Analytical solution
uana=A*c;                                                                           %Analytical solution
u0=uana(2);
[s,c12,d] = ellipj(sqrt(omeg0^2+gamma*A^2)*(k/2),(gamma*A^2)/(2*gamma*A^2+2*omeg0^2)); %Analytical solution at the time k*1/2
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

    klin(n-1)=(1/2)*((u0-um1)/k)^2;                 %Kinetic energy

    vlin(n-1)=(1/k^2)*(1-cos(omeg0*k))*um1*u0;  %Linear potential energy

    Hsav(n-1)=klin(n-1)+vlin(n-1)+(1/2)*(psim1)^2;  %Energy

    g0=(gamma*u0^3)/(sqrt(gamma*(u0^4)/2+eps));     %g vector (scalar on this case)
  
    up1=(1/(1+(g0^2*k^2)/4))*((2*cos(omeg0*k))*u0+(((g0^2*k^2)/4)-1)*um1-k^2*g0*psim1); %Exact linear solution + SAV
    
    psi0=psim1+(g0/2)*(-um1+up1);
    
    % Storing 
    g(n)=g0;
    psi(n)=psi0;
    u(n+1)=up1;

    % State swap
    psim1=psi0;
    um1=u0;
    u0=up1;
    
end

%% Plot of numerical and analytical solutions
figure
plot(t,u,LineWidth=3)
hold on
plot(t,uana,LineWidth=3,LineStyle="--")
legend("u","uana")
xlabel("Time(s)")
ylabel("Amplitude")
set(gca,'Fontsize',20)

%% Plot of the energy

figure
plot(th,(Hsav-Hsav(3))./Hsav,LineWidth=3)

