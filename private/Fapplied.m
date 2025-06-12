function [pext,mdShapes] = Fapplied(A,t,Ts,Nx,Ny,Phi,D,rho,Lz,Om,k,h,Nmodes,ExType,ForType)


Fext=zeros(1,Ts);
mdShapes=zeros(Ny+1,Nx+1,Nmodes);
pext=zeros(Ts,Nmodes);
T=max(t);

switch ForType

    case "periodic" % Periodic excitation
        env=(tanh(4*t)+tanh(4*((T-2)-t)))/2; % smoothing at the begining and end of the signal
        f0=Om(1)./(2*pi)-100;
        Fext=sin(2*pi*f0*t).*env(:)';

    case "chirp" %linear chirp between f0 and f1
        f0=Om(1)./(2*pi)-100;
        f1=Om(1)./(2*pi)+10;
        env=tanh(4*t); %smooth initial excitation
        Fext=chirp(t,f0,T,f1).*env(:)';

    case "impact" %raised cosine excitation
        T0=0.0005;
        Thw=0.0005;
        c0=D/(rho*Lz);
        for n=1:Ts
            if abs(n*k-T0) <Thw
                Fext(1,n)=c0*(1+cos(pi*(n*k-T0)/Thw));
            end
        end


end

switch ExType

    case "point" % Exictation on a single point point

        for m = 1 : Nmodes
            mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
            pext(:,m)=A*Fext*mdShapes(floor(end/4),floor(end/4),m);
        end

    case "distributed" % Distributed excitation on the plate

        for m = 1 : Nmodes
            mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
            pext(:,m)=A*Fext*trapzIntcalc(Phi(:,m),h,Nx,Ny);
        end


end