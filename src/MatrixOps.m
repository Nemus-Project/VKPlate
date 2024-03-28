%-----------------------------------------------------------------------
%                 Matrix Derivative Operator for VK plates
%                       Michele Ducceschi
%                       University of Bologna
%                           24 Mar 2024
%------------------------------------------------------------------------



clear all
%close all
clc

%------------------------------------------------------------------------
% custom plate parameters
Lx = 1 ;
Ly = 0.7 ;
nums=linspace(0.01,1,100);

hsp  = nums*sqrt(Lx*Ly)*0.1 ;
%------------------------------------------------------------------------
for i=1:100
    %tic
h=hsp(i);
%------------------------------------------------------------------------
% derived parameters
Nx = floor(Lx/h) ; 
Ny = floor(Ly/h) ; 
x = (0:Nx)*h ;
y = (0:Ny)*h ;
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Dx
Dx = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 2 : Nx

    Dx((m-1)*(Ny+1)+1:m*(Ny+1),(m-2)*(Ny+1)+1:(m-1)*(Ny+1)) = speye(Ny+1)*(-1/2/h) ;
    Dx((m-1)*(Ny+1)+1:m*(Ny+1),m*(Ny+1)+1:(m+1)*(Ny+1)) = speye(Ny+1)*(1/2/h) ;

end

Blk00 = speye(Ny+1,Ny+1)*(-3/2/h) ; 
Blk01 = speye(Ny+1,Ny+1)*(2/h) ; 
Blk02 = speye(Ny+1,Ny+1)*(-1/2/h) ; 

Dx(1:Ny+1,1:Ny+1) = Blk00 ;
Dx(1:Ny+1,Ny+2:2*Ny+2) = Blk01 ;
Dx(1:Ny+1,2*Ny+3:3*Ny+3) = Blk02 ;

Dx(end-Ny:end,end-Ny:end) = -Blk00 ;
Dx(end-Ny:end,end-2*Ny-1:end-Ny-1) = -Blk01 ;
Dx(end-Ny:end,end-3*Ny-2:end-2*Ny-2) = -Blk02 ;
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Dy
DyBlk = sparse(Ny+1,Ny+1) ;


for n = 2 : Ny 

    DyBlk(n,n-1) = -1/2/h ; DyBlk(n,n+1) = 1/2/h  ;

end

 DyBlk(1,1) = -3/2/h ; DyBlk(1,2) = 2/h ; DyBlk(1,3) = -1/2/h ;
 DyBlk(end,end) = 3/2/h ; DyBlk(end,end-1) = -2/h ; DyBlk(end,end-2) = 1/2/h ;

Dy = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 1 : Nx + 1

    Dy((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1)) = DyBlk ;
    
end
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% test out against some trig functions
[X,Y]   = meshgrid(x,y) ;

f       = cos(4*X/Lx).*sin(5*Y/Ly) ;
fx      = -4/Lx*sin(4*X/Lx).*sin(5*Y/Ly) ;
fy      = 5/Ly*cos(4*X/Lx).*cos(5*Y/Ly) ;
fxy     = -4/Lx*5/Ly*sin(4*X/Lx).*cos(5*Y/Ly) ;
fxx     = -4/Lx*4/Lx*cos(4*X/Lx).*sin(5*Y/Ly) ;
fyy     = -5/Ly*5/Ly*cos(4*X/Lx).*sin(5*Y/Ly) ;

%-- compute higher derivatives using composition. Note that Dx and Dy commute (as expected)
F       = reshape(f,[(Nx+1)*(Ny+1) 1]) ;
Fx      = Dx*F ;
Fy      = Dy*F ;
Fxy     = Dx*Dy*F ; %-- can use Fxy = Dy*Dx*F ; as well (commutative operation)
Fxx     = Dx*Dx*F ;
Fyy     = Dy*Dy*F ;
Dxx     = Dx*Dx   ;
Dyy     = Dy*Dy   ;
Dxy     = Dx*Dy   ;


% figure
% subplot(1,2,1)
% mesh(X,Y,fx) ; view(2); title('fx')
% subplot(1,2,2)
 Fx = reshape(Fx,[(Ny+1) (Nx+1)]) ;
% mesh(X,Y,Fx) ; view(2); title('Fx')
% 
% 
% figure
% subplot(1,2,1)
% mesh(X,Y,fy) ; view(2); title('fy')
% subplot(1,2,2)
 Fy = reshape(Fy,[(Ny+1) (Nx+1)]) ;
% mesh(X,Y,Fy) ; view(2); title('Fy')
% 
% 
% figure
% subplot(1,2,1)
% mesh(X,Y,fxy) ; view(2); title('fxy')
% subplot(1,2,2)
 Fxy = reshape(Fxy,(Ny+1),(Nx+1)) ;
% mesh(X,Y,Fxy) ; view(2); title('Fxy')
% 
% figure
% subplot(1,2,1)
% mesh(X,Y,fxx) ; view(2); title('fxx')
% subplot(1,2,2)
 Fxx = reshape(Fxx,(Ny+1),(Nx+1)) ;
% mesh(X,Y,Fxx) ; view(2); title('Fxx')
% 
% figure
% subplot(1,2,1)
% mesh(X,Y,fyy) ; view(2); title('fyy')
% subplot(1,2,2)
 Fyy = reshape(Fyy,(Ny+1),(Nx+1)) ;
% mesh(X,Y,Fyy) ; view(2); title('Fyy')
%------------------------------------------------------------------------
errmidx(i)=abs(fx(ceil(Ny/2),ceil(Nx/2))-Fx(ceil(Ny/2),ceil(Nx/2)));
errmidy(i)=abs(fy(ceil(Ny/2),ceil(Nx/2))-Fy(ceil(Ny/2),ceil(Nx/2)));
errmidxx(i)=abs(fxx(ceil(Ny/2),ceil(Nx/2))-Fxx(ceil(Ny/2),ceil(Nx/2)));
errmidxy(i)=abs(fxy(ceil(Ny/2),ceil(Nx/2))-Fxy(ceil(Ny/2),ceil(Nx/2)));
errmidyy(i)=abs(fyy(ceil(Ny/2),ceil(Nx/2))-Fyy(ceil(Ny/2),ceil(Nx/2)));

errbcx(i)=abs(fx(end,end)-Fx(end,end));
errbcy(i)=abs(fy(end,end)-Fy(end,end));
errbcxx(i)=abs(fxx(end,end)-Fxx(end,end));
errbcxy(i)=abs(fxy(end,end)-Fxy(end,end));
errbcyy(i)=abs(fyy(end,end)-Fyy(end,end));
%toc
end
%%
 close all
% 
% figure
% loglog(hsp,errmidxx,LineWidth=2)
% hold on
% loglog(hsp,errmidyy,LineWidth=2)
% loglog(hsp,errmidxy,LineWidth=2)
% loglog(hsp,errmidx,LineWidth=2)
% loglog(hsp,errmidy,LineWidth=2)
% legend("Dxx","Dyy","Dxy","Dx","Dy")
% xlabel("h")
% ylabel("Absolute error")
% set(gca,'Fontsize',20)
% title("Error in the middle",'FontSize',25)

figure
loglog(hsp,errbcxx,LineWidth=2)
hold on
loglog(hsp,errbcyy,LineWidth=2)
loglog(hsp,errbcxy,LineWidth=2)
loglog(hsp,errbcx,LineWidth=2)
loglog(hsp,errbcy,LineWidth=2)
legend("Dxx","Dyy","Dxy","Dx","Dy")
%legend("Dxx","Dyy","Dx","Dy")
xlabel("h")
ylabel("Absolute error")
set(gca,'Fontsize',20)
title("Error at the boundary",'FontSize',25)

%%
slope_i2_i1 = log10(errbcyy(1)/errbcyy(2))/log10(hsp(1)/hsp(2))

