clear all
close all
clc
%%
load("./param/modesVK.mat");
for m = 1 : Nmodes
    mdShapes(:,:,m) = reshape(V(:,m),[(Ny+1),(Nx+1)]) ;
   %subplot(4,4,m)
   mesh(3000*(mdShapes(:,:,m)),(abs(mdShapes(:,:,m))),'FaceColor','texturemap') ;  %axis equal; axis tight ; view(2) ;
end
%%
Id=speye(Nx+1,Ny+1);
Dxx=sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
ef = ones((Nx+1)*(Ny+1),1)/h^2;
Dxx=spdiags([ef -2*ef ef],[-Nx-1,0,Nx+1],(Nx+1)*(Ny+1),(Nx+1)*(Ny+1))

figure
spy(Dxx)
%clear Dxx
%Dxx=Dxx/(h^2);


test=Dxx*V(:,1);
for m = 1 : Nmodes
    testm = reshape(test(:,1),[(Ny+1),(Nx+1)]) ;
end
figure
mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;

%%
e = ones(Ny+1,1)/h^2;
dy=spdiags([e -2*e e],-1:1,Nx+1,Ny+1);
%dy(1,1)=3;
%dy(end,end)=3;
figure
spy(dy)
%%
Dyy=sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
Dyy(1:Nx+1,1:Ny+1)=dy;

for i =1:Nx
    Dyy(i*(Nx+1)+1:(i+1)*(Nx+1),i*(Ny+1)+1:(i+1)*(Ny+1))=dy;
end
%Dxx((Nx-1)*Nx+1:Nx*Nx,(Nx-1)*Ny+1:Nx*Ny)=-2*Id;
%Dxx((Nx+1)*Nx+1:(Nx+2)*Nx,(Nx+1)*Ny+1:(Nx+2)*Ny)=-2*Id;

figure
spy(Dyy)

test=Dyy*V(:,1);
for m = 1 : Nmodes
    testm = reshape(test(:,1),[(Ny+1),(Nx+1)]) ;
end
figure
mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;
%%


%%

Testvec=[1:(Nx+1)*(Ny+1)]';

for m = 1 : Nmodes
    testm = reshape(Testvec,[(Ny+1),(Nx+1)]) ;
end
figure
mesh(3000*(testm(:,:,1)),(abs(testm(:,:,1))),'FaceColor','texturemap') ;

%%
Dxy=fidimat(Ny+1,Nx+1,'grad',2);




