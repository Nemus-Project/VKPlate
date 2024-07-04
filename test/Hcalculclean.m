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

test=Dxx*V(:,nmod);
for m = 1 : Nmodes
    testm = reshape(test(:,1),[(Ny+1),(Nx+1)]) ;
end

figure
mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;

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
%%
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
maxana=max(abs(V(:,nmod)));
Ana=zeros((Ny+1),(Nx+1));

for i=1:length(x)
    for j=1 :length(y)
        Ana(i,j)=maxana*sin((nmod*pi/Lx)*x(i))*sin((mmod*pi/Ly)*y(j));
        hanam(i,j)=Ana(i,j)*((h^2*2*nmod^4*pi^4)/(Lx^2*Ly^2))*((Ana(i,j)^2)-(h^2*maxana*cos((nmod*pi/Lx)*x(i))*cos((mmod*pi/Ly)*y(j))^2));
    end
end

figure
mesh(3000*(Ana(:,:)),(abs(Ana(:,:))),'FaceColor','texturemap') ;
ana = reshape(Ana(:,:),[(Ny+1)*(Nx+1),1]) ;
%hanam=reshape(hana(:,1),[(Ny+1),(Nx+1)]) ;
%%
W=V;
Z=abs(V);

hv=Z(:,lmod).*((Dxx*V(:,nmod)).*(Dyy*W(:,mmod))+(Dxx*W(:,mmod)).*(Dyy*V(:,nmod))-2*(Dxy*V(:,nmod)).*(Dxy*W(:,mmod)));
hv=V(:,lmod).*((Dxx*V(:,nmod)).*(Dyy*V(:,mmod))+(Dxx*V(:,mmod)).*(Dyy*V(:,nmod))-2*(Dxy*V(:,nmod)).*(Dxy*V(:,mmod)));

%hana=ana(:,lmod).*((Dxx*ana(:,nmod)).*(Dyy*ana(:,mmod))+(Dxx*ana(:,mmod)).*(Dyy*ana(:,nmod))-2*(Dxy*ana(:,nmod)).*(Dxy*ana(:,mmod)));

for m = 1 : Nmodes
    hm = reshape(hv(:,1),[(Ny+1),(Nx+1)]) ;
    %hanam=reshape(hana(:,1),[(Ny+1),(Nx+1)]) ;
end
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
figurepyp
mesh(3000*(hm(:,:,1)),(abs(hm(:,:,1))),'FaceColor','texturemap') ;
figure
mesh(3000*(hanam(:,:,1)),(abs(hanam(:,:,1))),'FaceColor','texturemap') ;

%%

%%
X=linspace(0,Lx*Ly,(Nx+1)*(Ny+1));
H=trapz(y,trapz(x,hm,2))
Hana=trapz(y,trapz(x,hanam,2))

HH=h^2*Intmat*hv
HHana=h^2*Intmat*hana