clear all
close all
clc
%%
nbpts=150;
blkt=zeros(nbpts,1);
diagt=zeros(nbpts,1);
nxny=zeros(nbpts,1);
for n=nbpts:-1:1
Nx=10*n;
Ny=1000;
h=0.01;
nxny(n)= Nx*Ny;
tic
Dyy=DyyBuild(Nx,Ny,h,'blk');
blkt(n)=toc;

tic
Dyy=DyyBuild(Nx,Ny,h,'diag');
diagt(n)=toc;
disp(n)
end
%%
figure
plot(nxny,diagt,nxny,blkt,Linewidth=3)
legend("Diag","Block ")
xlabel("Number of operator points")
ylabel("Time (s)")
set(gca,'Fontsize',20)