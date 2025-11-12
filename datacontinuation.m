clear all
close all
clc

%% Importing modal parameters

load("./param/Test20modespurple50PsiT=1.mat")


load("./w10step2dB20modespurpleT=1w0=1e-9f=f0+20conti.mat")

%% Importing the time signal



Id=speye(Nmodes);

c0=0.08;                                            %Frequency independant damping term

c1=0.04;                                              % Frequency dependant damping term

xi=Id;                                     % Modal damping matrix

xi=(c1*Om.^(0.75)+c0).*Id;



%% extracting 1 period
Fs=44100;

figure
plot(wp)

time1=1192640;
time2=1192937;

Pextr=wp(time1:time2);%extract period
Pextr(end)=Pextr(1);%ensure a true period
P=(time2-time1)/Fs;%period duration

PextQ=Qsv(time1:time2,:);
PextEta=ETASV(time1:time2,:);
figure
plot(Pextr)
%% 
q=zeros(Nmodes,1);

a=zeros(Nmodes,1);
a(1)=1e-9;

eta=zeros(Npsi,1);

Hv4=Hv;
Hv2=Hv;


for p = 1:Nmodes
    for q = p:Nmodes
        for i = 1:Npsi
            Hv4(i,p,q)= Hv(i,p,q)/zetafourth(i);
            Hv4(i,q,p)= Hv(i,q,p)/zetafourth(i);
        end
    end
end


for p = 1:Nmodes
    for q = p:Nmodes
        for i = 1:Npsi
            Hv2(i,p,q)= Hv(i,p,q)/sqrt(zetafourth(i));
            Hv2(i,q,p)= Hv(i,q,p)/sqrt(zetafourth(i));
        end
    end
end

Hv0= reshape(Hv,[Npsi*Nmodes,Nmodes]);
Hv2 = reshape(Hv2,[Npsi*Nmodes,Nmodes]); % reshape of the tensor into a matrix
Hv4 = reshape(Hv4,[Npsi*Nmodes,Nmodes]); % reshape of the tensor into a matrix


%% save the file
save("purplecontf0+20Hz2ndregime.mat","xi","Hv","Hv0","Hv4","Om","PextEta","PextQ","Pextr","Qsv","ETASV","a")


%% compute the coupling term from eta and q (and a)
 t5 = Hv0*(q+a);
 t5 = reshape(t5,[Npsi,Nmodes]);
 Gb = (t5.'*eta);

%% compute eta from q

 % Hks/zeta^4 * qi*qj

    t7 = Hv4*q;

    t7 = reshape(t7,[Npsi,Nmodes]);
    Gc_1 = t7*q;

    % Hks/zeta^4 * a * 2q
    t8 = Hv4*a;
    t8 = reshape(t8,[Npsi,Nmodes]);
    Gc_2 = t8*2*q;

    Gc = (Gc_1 + Gc_2);

    etasvs =- (E*Lz)*Gc;
%%