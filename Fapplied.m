function [pext,mdShapes] = Fapplied(A,t,Ts,Nx,Ny,Phi,D,rho,Lz,Om,k,h,Nmodes,ExType,ForType,f0)


Fext=zeros(1,Ts);
mdShapes=zeros(Ny+1,Nx+1,Nmodes);
pext=zeros(Ts,Nmodes);
T=max(t);

switch ForType

    case "periodic" % Periodic excitation
        env=(tanh(4*t)+tanh(4*((T-1)-t)))/2; % smoothing at the begining and end of the signal
        %env=1;
        %f0=Om(1)./(2*pi);
        Fext=sin(2*pi*f0*t).*env(:)';

    case "chirp" %linear chirp between f0 and f1
        %f0=Om(1)./(2*pi)-100;
        %f1=Om(1)./(2*pi)+100;
        f0=100;
        f1=2500;
        env=tanh(4*t); %smooth initial excitation
        Fext=chirp(t,f0,T,f1).*env(:)';

    case "chirpreverse" %linear chirp between f0 and f1
        %f0=Om(1)./(2*pi)-100;
        %f1=Om(1)./(2*pi)+10;
        f1=700;
        env=tanh(4*t); %smooth initial excitation
        Fext=fliplr(chirp(t,f0,T,f1).*env(:)');

    case "impact" %raised cosine excitation
        T0=0.0005;
        Thw=0.0005;
        c0=D/(rho*Lz);
        for n=1:Ts
            if abs(n*k-T0) <Thw
                Fext(1,n)=c0*(1+cos(pi*(n*k-T0)/Thw));
            end
        end

    case "ampstep"

        nbstep=3;% number of step in the signal
        per=T/nbstep;% Time of one step
        k=10 ; % Parameter that controls smoothness of transition,
        dl=0.2; % Parameter that delay the nonzero signal to be smoother
        %ex= [0.06 0.2 0.6 2.0 6.0 10];
        vini=0.1; %First voltage
        ex=vini*ones(nbstep,1);
        for loo=2:1:nbstep
            ex(loo)=ex(loo-1)*sqrt(2); %sqrt(2) corresponds to a 3dB amplitude step in excitation
        end


        st=0; %Initialisation of the enveloppe
        for exloo=1:1:nbstep
            if exloo>1
                st=st+(ex(exloo)-ex(exloo-1))*((1+tanh(k*(t-(exloo-1)*per-dl)))/2);
            else
                st=st+ex(exloo)*((1+tanh(k*(t-(exloo-1)*per-dl)))/2);
            end
        end

        %f0=Om(1)./(2*pi);
        Fext=st.*sin((2*pi*f0)*t);




end

switch ExType

    case "point" % Exictation on a single point point

        for m = 1 : Nmodes
            mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
            pext(:,m)=A*Fext*mdShapes(floor(end/2),floor(end/2),m);
        end

    case "distributed" % Distributed excitation on the plate

        for m = 1 : Nmodes
            mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
            pext(:,m)=A*Fext*trapzIntcalc(Phi(:,m),h,Nx,Ny);
        end

    case "tilted"

        for m = 1 : Nmodes
            tilt=linspace(-0.002,0.002,(Ny+1));
            tilt2D=ones((Ny+1),(Nx+1)).*tilt;
            tiltsolve=reshape(tilt2D,[1,(Ny+1)*(Nx+1)])';
            mdShapes(:,:,m) = reshape(Phi(:,m),[(Ny+1),(Nx+1)]) ;
            pext(:,m)=A*Fext*trapzIntcalc(tiltsolve.*Phi(:,m),h,Nx,Ny);
        end
end