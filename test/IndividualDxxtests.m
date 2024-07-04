clear all
close all
clc
%%
load("./param/modesint.mat");
Dxy=fidimat(Ny+1,Nx+1,'grad',1);
Dxx=fidimat(Ny+1,Nx+1,'xx',1);
Dyy=fidimat(Ny+1,Nx+1,'yy',1);

Intmat=ones(1,(Nx+1)*(Ny+1));
Intmat([2:Nx,end-Nx+1:end-1])=1/2;
for i=Nx+2:((Nx+1)*(Ny+1))-Nx+1
    if mod(i,Nx+1)==0
        Intmat(i)=1/2;
    end
    if mod(i,Nx+1)==1
        Intmat(i)=1/2;
    end

end
Intmat([1,Nx+1,end-Nx,end])=1/4;
test2=reshape(Intmat,[(Ny+1),(Nx+1)]) ;

nmod=1;
mmod=1;
lmod=1;



% 
% test=Dyy*V(:,nmod);
% for m = 1 : Nmodes
%     testm = reshape(test(:,1),[(Ny+1),(Nx+1)]) ;
% end
% figure
% mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;
% 
% test=Dxy*V(:,nmod);
% for m = 1 : Nmodes
%     testm = reshape(test(:,1),[(Ny+1),(Nx+1)]) ;
% end
% figure
% mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;
%%
for m = 1 : Nmodes
    mode1 = reshape(V(:,1),[(Ny+1),(Nx+1)]) ;
end
figure
mesh(3000*(mode1(:,:)),(abs(mode1(:,:))),'FaceColor','texturemap') ;
title "mode1"
%%
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
maxana=(-min(V(:,nmod)) < max(V(:,nmod)))*(max(V(:,nmod))-min(V(:,nmod)))+min(V(:,nmod));
Ana=zeros((Ny+1),(Nx+1));

for i=1:length(x)
    for j=1 :length(y)
        Ana(i,j)=maxana*sin((nmod*pi/Lx)*x(i))*sin((mmod*pi/Ly)*y(j));
        Deriv(i,j)=maxana*(nmod^2*pi^2/(Ly*Lx))*cos((nmod*pi/Lx)*x(i))*cos((mmod*pi/Ly)*y(j));
        hanam(i,j)=Ana(i,j)*((2*nmod^4*pi^4)/(Lx^2*Ly^2))*((Ana(i,j)^2)-(maxana*cos((nmod*pi/Lx)*x(i))*cos((mmod*pi/Ly)*y(j))^2));
        
    end
end
Deriv([1,end],:)=0;
Deriv(:,[1,end])=0;
figure
mesh(3000*(Ana(:,:)),(abs(Ana(:,:))),'FaceColor','texturemap') ;
title "mode1ana"
ana = reshape(Ana(:,:),[(Ny+1)*(Nx+1),1]) ;
deriv = reshape(Deriv(:,:),[(Ny+1)*(Nx+1),1]) ;
%hanam=reshape(hana(:,1),[(Ny+1),(Nx+1)]) ;
%%

test=Dxy*V(:,nmod);
for m = 1 : Nmodes
    testm = reshape(test(:,1),[(Ny+1),(Nx+1)]) ;
end

figure
mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;
title "DxyV"

%testana=Dxx*ana;
%testana=-(h^2*nmod^2*pi^2/(Ly*Lx))*ana;
testana=deriv;
testanam = reshape(testana(:,1),[(Ny+1),(Nx+1)]) ;

figure
mesh(3000*(testanam(:,:,1)),(abs(testanam(:,:,1))),'FaceColor','texturemap') ;
title "Dxyana"

testdif=(testana-test);

ermax=max(abs(testdif))
