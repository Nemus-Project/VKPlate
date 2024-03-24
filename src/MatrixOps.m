%-----------------------------------------------------------------------
%                 Matrix Derivative Operator for VK plates
%                       Michele Ducceschi
%                       University of Bologna
%                           24 Mar 2024
%------------------------------------------------------------------------



clear all
close all
clc

%------------------------------------------------------------------------
% custom plate parameters
Lx = 1 ;
Ly = 0.7 ;
h  = sqrt(Lx*Ly)*0.01 ;
%------------------------------------------------------------------------


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

figure
subplot(1,2,1)
mesh(X,Y,fx) ; view(2); title('fx')
subplot(1,2,2)
Fx = reshape(Fx,[(Ny+1) (Nx+1)]) ;
mesh(X,Y,Fx) ; view(2); title('Fx')


figure
subplot(1,2,1)
mesh(X,Y,fy) ; view(2); title('fy')
subplot(1,2,2)
Fy = reshape(Fy,[(Ny+1) (Nx+1)]) ;
mesh(X,Y,Fy) ; view(2); title('Fy')


figure
subplot(1,2,1)
mesh(X,Y,fxy) ; view(2); title('fxy')
subplot(1,2,2)
Fxy = reshape(Fxy,(Ny+1),(Nx+1)) ;
mesh(X,Y,Fxy) ; view(2); title('Fxy')

figure
subplot(1,2,1)
mesh(X,Y,fxx) ; view(2); title('fxx')
subplot(1,2,2)
Fxx = reshape(Fxx,(Ny+1),(Nx+1)) ;
mesh(X,Y,Fxx) ; view(2); title('Fxx')

figure
subplot(1,2,1)
mesh(X,Y,fyy) ; view(2); title('fyy')
subplot(1,2,2)
Fyy = reshape(Fyy,(Ny+1),(Nx+1)) ;
mesh(X,Y,Fyy) ; view(2); title('Fyy')
%------------------------------------------------------------------------










