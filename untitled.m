clear all
close all
clc 
%%

load("./param/Test6modesgreen.mat")
Ev=permute(Hv,[3 2 1]);
Gamtest=zeros(6,6,6,6);%for loop test
%Gam=zeros(6,6,6,6);% tensor product test

%% for loop
% for s= Nmodes
%     for k = 1: Nmodes
%         for  m= 1: Nmodes
%             for n = 1: Nmodes
%                 for l=1:Npsi
%                     Gamtest(s,k,m,n)=Gamtest(s,k,m,n)+ (Hv(l,m,n)*Ev(s,k,l))/2*zetafourth(l);
%                 end
%             end
%         end
%     end
% end

for s= 1:Nmodes
    for k = 1: Nmodes
        for  m= 1: Nmodes
            for n = 1: Nmodes
                for l=1:Npsi
                    Gamtest(s,k,m,n)=Gamtest(s,k,m,n)+ (Hv(l,k,m)*Hv(l,n,s))/2*zetafourth(l);
                end
            end
        end
    end
end

%% tensorproduct
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
% 
% for p = 1:Nmodes
%     for q = p:Nmodes
%         for i = 1:Npsi
%             Hvtest(i,p,q)= Hv4(i,p,q)*zetafourth(i);
%             Hvtest(i,q,p)= Hv4(i,q,p)*zetafourth(i);
%         end
%     end
% end



Gam = tensorprod(Hv,Hvmod,1,1);

%%

idtest=Gam-Gamtest;

%% test multiplications
q0=zeros(Nmodes,1);
a0=zeros(Nmodes,1);
q0=[1:Nmodes]'; %vertical vector
qs=[1:Nmodes]'; %

eta=[1:Npsi]'; %vertical vector
etas=[1:Npsi]';
etam=[1:Npsi]';

Hv0= reshape(Hv,[Npsi*Nmodes,Nmodes]);
Hv2 = reshape(Hv2,[Npsi*Nmodes,Nmodes]); % reshape of the tensor into a matrix
Hv4 = reshape(Hv4,[Npsi*Nmodes,Nmodes]); % reshape of the tensor into a matrix
C=eye(6);


%left hand term
t0=Hv2*(q0+a0);%half of the sum
t0 = reshape(t0,[Npsi,Nmodes]);%reshape into a matrix
G0 = (t0.'*t0);%Matrix of the left hand double sum


mat_imp = C + G0; % implicit matrix


%right hand term

 t1 = Hv2*q0;%half of the sum
 t1 = reshape(t1,[Npsi,Nmodes]);
 t4 = t1*a0;
 Ga = t0.'*t4; % Matrix

eta_temp = (etam-eta);
    t5 = Hv0*(q0 + a0);
    t5 = reshape(t5,[Npsi,Nmodes]);
    Gb = (t5.'*eta_temp);
    
    Gtotal = Ga - Gb;





%Eta vector

t7 = Hv4*q0;
t7 = reshape(t7,[Npsi,Nmodes]);
Gc_1 = t7*q0;


 t8 = Hv4*a0;
 t8 = reshape(t8,[Npsi,Nmodes]);
 Gc_2 = t8*(qs + q0);

Gc = (Gc_1 + Gc_2);
    
etas = -eta - Gc;
%Update
 qm = q0;
 q0 = qs;
 etam = eta;
 eta = etas;
