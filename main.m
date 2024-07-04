clear all
close all
clc

% This is the main code of VKPlate. It calls a parameter file generated
% by genparams.m, and use an energy exact scheme to
% compute all the modal displacement of a nonlinear plate. It also has a
% comparison with Stormer-Verlet integration scheme.
%
% The input file should include :
%
%   USEFUL FOR CONTROL
%       BCsPhi - Displacement boundary conditions   | 4 x 2 vector
%       BCsPsi - Stress boundary conditions         | 4 x 2 vector
%       Om2 - Angular frequency of the Psi modes    | Nmodes x 1 vector
%       Lx - x length of the plate                  | float
%       Ly - y length of the plate                  | float
%       X - x space                                 | Ny+1 x Nx+1 matrix
%       Y - y space                                 | Ny+1 x Nx+1 matrix
%
%   USED FOR COMPUTATION
%       E - Young modulus                           | float
%       rho - Density                               | float
%       nu - Poisson's ratio                        | float
%       h - Space step                              | float
%       Lz - Space step                             | float
%       Nmodes  - Number of modes                   | float
%       Hv - Coupling coeficients                   | Nmodes x Nmodes x Nmodes matrix
%       Nx - Number of points in the x direction    | float
%       Ny - Number of points in the y direction    | float
%       Om - Angular frequency of the Phi modes     | Nmodes x 1 vector
%       Phi - Eigenvectors of displacement          | (Nx+1*Ny+1) x Nmodes matrix
%       Psi - Eigenvectors of stress                | (Nx+1*Ny+1) x Nmodes matrix
%       zetafourth - Eigenvalues of the stress      | Nmodes x 1 vector
%
% The raw output is the Q matrix, that can be transformed by to physical
% displacement (W) by projecting it on the Phi vector.


%% Loading parameter file

load("./param/Test_NL_Fullclamp_3.mat"); %Parameter file is generated by main

%% Declared parameters

Fs      = 44100 ;                     % Sampling rate
os      = 2 ;                         % Oversampling factor
T       = 2;                          % Total time

%% Derived variable

D       = E * Lz^3 / 12 / (1-nu^2) ;
k       = (1/Fs)/os ;               % Timestep
Ts      = floor(T/k) ;              % Number of time grid points
t       =linspace(0,T,Ts);          % Time vector
x       =linspace(0,Lx,Nx+1);       % x space vector
y       =linspace(0,Ly,Ny+1);       % y space vector


%% Defining excitiation force

Fext=zeros(1,Ts);
c0=D/(rho*Lz);
T0=0.0005;
Thw=0.0005;
for n=1:Ts
    if abs(n*k-T0) <Thw
        Fext(1,n)=0.5*c0*(1+cos(pi*(n*k-T0)/Thw));
    end
end


for m = 1 : Nmodes
    mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
    pext(:,m)=0e-1*Fext*mdShapes(floor(end/4),floor(end/4),m);

end
pext=pext';

%% Variable declaration

Id=speye(Nmodes);                               % Identity matrix

Dw1=diag((Om(1:Nmodes)));                       % Diagonal matrix of Omega

Dw2=diag((Om(1:Nmodes)).^2);                    % Diagonal matrix of Omega^2

Dw=diag((cos(k*Om(1:Nmodes))));                 % Diagonal for the exact lossless integrator

DV0=2*k^(-2)*(Id-diag((cos(k*Om(1:Nmodes)))));  % 

Ddiag=2*sparse(Dw);                             % 2*Sparse of Dw

D0sv=2*Id-(k^2)*Dw2;                            % 


qm=zeros(Nmodes,1);                             % Modal displacement vector at n-1

q0=zeros(Nmodes,1);                             % Modal displacement vector at n

qmsv=zeros(Nmodes,1);                           % Modal displacement vector at n-1 Stormer Verlet

q0sv=zeros(Nmodes,1);                           % Modal displacement vector at n Stormer Verlet

qssv=zeros(Nmodes,1);                           % Modal displacement vector at n+1 Stormer Verlet

nlsv=zeros(Nmodes,1);                           % Nonlinear term of Stormer Verlet

dUdq=zeros(Nmodes,1);                           % Derivative of the energy on q

g=zeros(Nmodes,1);                              % SAV nonlinear vector

En=zeros(1,Ts);                                 % Kinetic Energy

V0=zeros(1,Ts);                                 % Linear Potential

PSI0=zeros(1,Ts);                               % Psi time vector

Q=zeros(Ts,Nmodes);                             % q time vector

Qsv=zeros(Ts,Nmodes);                           % q time vector Stormer-Verlet

ETA=zeros(Ts,Nmodes);                           % eta time vector

eta=zeros(Nmodes,1);                            % Modal stress vector at n

etam=zeros(Nmodes,1);                           % Modal stress vector at n-1

xi=1e-1*Id;                                     % Damping matrix

DSM=(0.5*k*xi*Dw1+Id)^-1;                       % Damping term in the scheme

Dqmsv=0.5*k*xi*Dw1-Id;                          % Damping term in the Stormer-Verlet scheme


%% Initialisation

%qm=flipud(sort(1e-4*rand(Nmodes,1)));           %random excitation

qm(1)=5e-3;                                     % Initialisation of the displacement at n-1

q0=qm;                                          % Initialisation of the displacement at n

qmsv=qm;                                        % Initialisation of the displacement at n-1 Stormer-Verlet

q0sv=qmsv;                                      % Initialisation of the displacement at n Stormer-Verlet

Q(1,:)=q0;

Qsv(1,:)=q0sv;

enerr=1/(rho*Lz);                               % Energy multiplicator (see paper)

for ki = 1:Nmodes
    for ji =1:Nmodes
        for ii =1:Nmodes
            eta(ki)=eta(ki)-((E*Lz)/(2*zetafourth(ki)))*Hv(ii,ji,ki)*q0(ii)*q0(ji); % Initialisation of eta
        end
    end
end

etasv=eta;

ETA(1,:)=eta;

ETASV(1,:)=etasv;

U=enerr*(1/(2*E*Lz))*zetafourth(1:Nmodes)'*0.5*(eta.^2+eta.^2);   %to change with actual energy (meaning etan-1)

ep=1e-10;                                       % Potential shift

psim1=sqrt(2*U+ep);                             % Initialisation of psi at n-1/2

%% Integration loop

tic
for n = 2 : Ts

    En(n)=0.5*sum(((q0-qm)/k).^2);  %Kinetic energy at timestep n

    V0(n)=0.5*sum(DV0*q0.*qm);      %Linear potential energy at timestep n

    eta(:)=0;
    etasv(:)=0;
    for li = 1:Nmodes
        for mi =1:Nmodes
            for ni =1:Nmodes
                eta(li)=eta(li)-((E*Lz)/(2*zetafourth(li)))*Hv(mi,ni,li)*q0(mi)*q0(ni);         %eta computation 
                etasv(li)=etasv(li)-((E*Lz)/(2*zetafourth(li)))*Hv(mi,ni,li)*q0sv(mi)*q0sv(ni); %eta computation for Stormer Verlet

            end
        end
    end

    ETA(n,:)=eta;

    nlsv(:)=0;
    for ki = 1:Nmodes
        for ji =1:Nmodes
            for ii =1:Nmodes
                nlsv(ki)=nlsv(ki)+Hv(ki,ji,ii)*q0sv(ii)*etasv(ji);      %Stormer Verlet nonlinear term computation
            end
        end
    end
 
    U=enerr*(1/(2*E*Lz))*zetafourth(1:Nmodes)'*(eta.^2);        %Nonlinear potential energy computation

    dUdq(:)=0;
    for ji = 1:Nmodes
        for ki =1:Nmodes
            for ii =1:Nmodes
                dUdq(ji)=dUdq(ji)-eta(ki)*(q0(ii)*Hv(ii,ji,ki));    %dU/dq computation

            end
        end
    end
    g=enerr*(1/sqrt(2*U+ep))*dUdq;  %g vector computation at each timestep

    G(n)=g(1);

    Q(n,:)=q0;

    Qsv(n,:)=q0sv;

    % Matrix definitions

    Dqnp=0.25*(k^2)*g*g'+Id;

    Dqnm=0.25*(k^2)*g*g'-Id;

    Dexp=Id-(0.25*(k^2)*g*g')/(1+0.25*k^2*g'*g);

    Dqnm2=0.25*(k^2)*g*g'+0.5*xi*Dw1*k-Id;
   
    Dexp2=DSM-(DSM*0.25*(k^2)*g*g'*DSM)/(1+0.25*k^2*g'*DSM*g);
    
    % Timestepping scheme

    qs=Dexp2*(D0sv*q0+Dqnm2*qm-g*(k^2)*psim1);% +k^2*pext(:,n)); %SAV

    qssv=DSM*(D0sv*q0sv+Dqmsv*qmsv+((k^2)/(rho*Lz))*nlsv); %Stormer Verlet

    psi0=psim1+0.5*g'*(qs-qm); % Psi time stepping




    % w(:,:)=0;
    % for m = 1 : Nmodes
    %     w=w+q0(m)*mdShapes(:,:,m);
    % end
    %W(:,:,n)=w;

    % State swap

    PSI0(n)=psim1;
    qm=q0;
    q0=qs;
    psim1=psi0;
    qmsv=q0sv;
    q0sv=qssv;


end
toc


%% Plotting the energy components
figure
plot(t,En-En(3),t,V0-V0(3),t,((PSI0.^2)-(PSI0(3)^2))/2)
figure
EN=En+V0+(PSI0.^2)/2;
scatter(t,(EN-EN(150))/EN(150),'.')
xlim([0.02 T])

%% Plotting Eta and Q
figure
plot(ETA(:,:),LineWidth=3)
title('\eta_k')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
figure
plot(Q(:,:),LineWidth=3)
title('w_k SAV')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
figure
plot(Qsv(:,:),LineWidth=3)
title('w_k Stormer-Verlet')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
%% Plotting the error
figure
plot(Q-Qsv,LineWidth=3)
title('w_k error')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)

%% Plotting psi
figure
plot(PSI0,LineWidth=3)
title('Psi')
xlabel('Timestep')
ylabel('Amplitude')
set(gca,'FontSize',20)
