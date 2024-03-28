function Dyy = DyyBuild(Nx,Ny,h,method)
% Dyy = DYYBUILD(Nx,Ny,h,method)
%
%       Nx     : number of grid points x-axis
%       Ny     : number of grid points y-axis
%       h      : grid spacing 
%       method : method for building matrix: 'blk' or 'diag'
%
%------------------------------------------------------------------------

% Parse args
if nargin < 4
    if (Nx*Ny) < 1000
        method = 'blk';
    else
        method = 'diag';
    end
end

%% Validate Arguments
valid_methods = ["blk", "diag"];

validateattributes(Nx, {'numeric'}, {'integer','positive'});
validateattributes(Ny, {'numeric'}, {'integer','positive'});
validateattributes(h,  {'numeric'}, {'real','positive'});
method = validatestring (method,    valid_methods);
%------------------------------------------------------------------------
% Dy
DyyBlk = sparse(Ny+1,Ny+1) ;


for n = 2 : Ny

    DyyBlk(n,n-1) = 1/h^2 ; DyyBlk(n,n+1) = 1/h^2  ; DyyBlk(n,n) = -2/h^2  ;

end

DyyBlk(1,1) = 2/h^2 ; DyyBlk(1,2) = -5/h^2 ; DyyBlk(1,3) = 4/h^2 ; DyyBlk(1,4) = -1/h^2 ;
DyyBlk(end,end) = 2/h^2 ; DyyBlk(end,end-1) = -5/h^2 ; DyyBlk(end,end-2) = 4/h^2 ; DyyBlk(end,end-3) = -1/h^2 ;

Dyy = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 1 : Nx + 1

    Dyy((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1)) = DyyBlk ;

end
%------------------------------------------------------------------------


end
