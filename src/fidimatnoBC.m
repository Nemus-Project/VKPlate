function outM = fidimatnoBC(l,m,h,type)
% FIDIMATNOBC Generate Finite Difference Spatial Sparcity Matrices
  %   FIDIMATNOBC(l,m,ord,bctype) A function to generate the biharmonic
  %   coefficients for 2D given the number of ALL grid points given by l and m.
  %
  %   bctype denotes boundary condition
  %
  %         Returns a Sparse matrix BH
  %
  %         
  %         l       % number of total grid points Y axis
  %         m       % number of total grid points X axis
  %         h       %
  %         type    % order of the matrix (string)
  %
  %                 % Valid type inputs
  %                 % 'Dxx', 'Dyy', 'Dxy', 'Dx', 'Dy'


%% derived parameters
Nx = l-1 ; 
Ny = m-1 ; 
x = (0:Nx)*h ;
y = (0:Ny)*h ;

%% Dx
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


%% Dy
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

switch type

    case 'Dx'
        outM = Dx;
    case 'Dy'
        outM = Dy;
    case'Dxx'
        outM = Dx*Dx;
    case 'Dyy' 
        outM = Dy*Dy;
    case 'Dxy'
        outM = Dx*Dy;
end
    
