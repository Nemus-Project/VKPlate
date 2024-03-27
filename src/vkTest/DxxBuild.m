function Dxx = DxxBuild(Nx,Ny,h)

%------------------------------------------------------------------------
% Dx
Dxx = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 2 : Nx

    Dxx((m-1)*(Ny+1)+1:m*(Ny+1),(m-2)*(Ny+1)+1:(m-1)*(Ny+1)) = speye(Ny+1)*(1/h^2) ;
    Dxx((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1))     = speye(Ny+1)*(-2/h^2) ;
    Dxx((m-1)*(Ny+1)+1:m*(Ny+1),m*(Ny+1)+1:(m+1)*(Ny+1))     = speye(Ny+1)*(1/h^2) ;

end

Blk00 = speye(Ny+1,Ny+1)*(2/h^2) ;
Blk01 = speye(Ny+1,Ny+1)*(-5/h^2) ;
Blk02 = speye(Ny+1,Ny+1)*(4/h^2) ;
Blk03 = speye(Ny+1,Ny+1)*(-1/h^2) ;

Dxx(1:Ny+1,1:Ny+1) = Blk00 ;
Dxx(1:Ny+1,Ny+2:2*Ny+2) = Blk01 ;
Dxx(1:Ny+1,2*Ny+3:3*Ny+3) = Blk02 ;
Dxx(1:Ny+1,3*Ny+4:4*Ny+4) = Blk03 ;

Dxx(end-Ny:end,end-Ny:end) = Blk00 ;
Dxx(end-Ny:end,end-2*Ny-1:end-Ny-1) = Blk01 ;
Dxx(end-Ny:end,end-3*Ny-2:end-2*Ny-2) = Blk02 ;
Dxx(end-Ny:end,end-4*Ny-3:end-3*Ny-3) = Blk03 ;
%------------------------------------------------------------------------

end
