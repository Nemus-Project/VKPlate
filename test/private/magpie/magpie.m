function [Om,Q,Nx,Ny,biHarm,Dm] = magpie(rho,E,nu,ldim,h,BCs,Nm,plot_type,shouldNormalise)
% MAGPIE What does this do?
%   [Om,Q,Nx,Ny,biHarm] = MAGPIE (density, Youngs, poisson, dim, h, BCs, Number_of_modes, plot_type)
%   A function that returns:
%           Om      : Angular modal frequencies
%           Q       : A matrix of column eigenvector(s)
%           Nx      : Grid points along the x-axis
%           Ny      : Grid points along the y-axis
%           biHarm  : Biharmonic Matrix for the plate
%           Dm      : the eigenvalues of the biharmonic matrix 
%
%   Arguments:
%       rho          %-- density [kg/m^3]
%       E            %-- Young's mod [Pa]
%       nu           %-- poisson's ratio
%
%       3 element array  representing  [x,y,z] dimensions of plate
%       ldim = [Lx,  %-- length along x [m]
%               Ly,  %-- length along y [m]
%               Lz]  %-- thickness [m]
%
%       h            %-- grid spacing
%
%       2 column array of eleastic boundary constants around each edge of
%       the plate.
%
%       Column 1 repesents  ...
%       Column 2 repesents  ...
%
%       BCs = [K0y, R0y;
%              Kx0, Rx0;
%              KLy, RLy;
%              KxL, RxL];
%
%       Nm      %-- number of modes to compute
%
%
%       Example:
%
%           %% physical and elastic parameters
%           Lx = 0.10; Ly = 0.08; Lz = 0.81e-3;
%           ldim = [Lx Ly Lz];   % plate dimensions [x, y, z] in metres
%           E       = 1.01e+11 ;        %-- Young's mod [Pa]
%           rho     = 8765 ;            %-- density [kg/m^3]
%           nu      = 0.3 ;             %-- poisson's ratio
%           Nm  = 16;               %-- number of modes to compute
%           h       = sqrt(Lx*Ly)*0.01; %--
%           BCs = ones(4,2) * 1e15      %-- elastic constants around the edges
%
%           [Om,Q,Nx,Ny,biHarm,Dm] = magpie(rho, E, nu, ldim, h, BCs, Nm,'none');
%% Varargs
if nargin < 9
    shouldNormalise = true;
end
if nargin < 8
    plot_type = 'none';
end
if nargin < 7
    Nm = 0;
end
%% Validation
validateattributes(rho,      {'double'}, {'nonempty'});
validateattributes(E,        {'double'}, {'nonempty'});
validateattributes(nu,       {'double'}, {'nonempty'});
validateattributes(ldim,     {'double'}, {'numel', 3});
validateattributes(h,        {'double'}, {'nonempty'});
validateattributes(BCs,      {'double'}, {'size', [4,2]});
validateattributes(Nm,   {'numeric'}, {'integer','nonnegative'});
validatestring(plot_type,["chladni","3D","none"]);

%% Unpack array variables
pack_ldim = num2cell(ldim);
pack_BCs = num2cell(BCs);
[Lx, Ly, Lz] = pack_ldim{:};
[K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};

%%--- derived parameters (don't change here)
D = E * Lz^3 / 12 / (1-nu^2);
% Nx      = floor(Lx/h) ;
% Ny      = floor(Ly/h) ;

Nx      = round(Lx/h) ;
Ny      = round(Ly/h) ;
%%----------------------------
%% Build BiHarmonic
biHarm = bhmat(BCs,[Nx Ny], h, Lz, E, nu);

%% EIGENVALUES

Nmodes = (Nx+1)*(Ny+1) ;
if Nm
    Nmodes = Nm ;
end
[Q,Dm] = eigs(biHarm,Nmodes,'smallestabs') ;
[~,indSort] = sort(diag((Dm))) ;
Q = Q(:,indSort) ;

Dm    = diag(Dm) ;
Om    = sqrt(abs(Dm))*sqrt(D/rho/Lz) ;
%freqs = Om/2/pi ;

if shouldNormalise
    for nQ = 1 : Nmodes
        Qtemp   = Q(:,nQ) ;
        Qnorm   = trapzIntcalc(Qtemp.*Qtemp,h,Nx,Ny) ;
        Qtemp   = Qtemp / sqrt(Qnorm) ;
        Q(:,nQ) = Qtemp ;
    end
end


switch plot_type
    case 'chladni'
        subs = ceil(sqrt(Nmodes));

        colormap('copper') ;
        cmp = colormap;
        cmp = flipud(cmp);
        colormap(cmp);

        for m = 1 : Nmodes
            mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
            subplot(subs,subs,m)
            mesh(3e3*real(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;
            view(2);
            axis equal;
            axis tight;
            axis off;
            clim([0.00005 0.002]);
        end

    case '3D'

        subs = ceil(sqrt(Nmodes));
        colormap('parula') ;
        xax = (0:Nx)*h ;
        yax = (0:Ny)*h ;
        [X,Y] = meshgrid(xax,yax) ;

        for m = 1 : Nmodes
            mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
            subplot(subs,subs,m)
            mesh(X,Y,3000*(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;
        end
end

end
