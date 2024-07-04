function [Qsign] = eigensign(Qref,Nxref,Nyref,Q,Nx,Ny,h,Nmodes,Lx,Ly)

xref=linspace(0,Lx,Nxref+1);
yref=linspace(0,Ly,Nyref+1);
[XREF,YREF]=meshgrid(xref,yref);
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
[X,Y]=meshgrid(x,y);
h=1;
Qsign=Q;
 for nQ = 1 : Nmodes
        Qtemp   = Q(:,nQ) ;

        Sptemp = reshape(Qtemp,[(Ny+1),(Nx+1)]) ;

        Qreftemp= Qref(:,nQ) ;

        Spreftemp=reshape(Qreftemp,[(Nyref+1),(Nxref+1)]) ;

        Spreftemp= interp2(XREF,YREF,Spreftemp,X,Y);
       

        Qinterp= reshape(Spreftemp,[(Ny+1)*(Nx+1),1]) ;
        % figure
        % subplot(2,2,1)
        % mesh(Spreftemp)
        % subplot(2,2,2)
        % mesh(Sptemp)
       
        Q1=Spreftemp-Sptemp;
        
        Q1vec= reshape(Q1,[(Ny+1)*(Nx+1)],1) ;

        Q2=Spreftemp+Sptemp;

        Q2vec= reshape(Q2,[(Ny+1)*(Nx+1)],1) ;

        % subplot(2,2,3)
        % mesh(Q1)
        % subplot(2,2,4)
        % mesh(Q2)
        %intq1=trapzIntcalc(Q1vec.*Q1vec,h,Nx,Ny);
        %intq2=trapzIntcalc(Q2vec.*Q2vec,h,Nx,Ny);
       
        sgn   = sign(-abs(trapzIntcalc(Q1vec.*Q1vec,h,Nx,Ny))+abs(trapzIntcalc(Q2vec.*Q2vec,h,Nx,Ny))); 
       
        Qtemp   = sgn*Qtemp;
        Qsign(:,nQ) = Qtemp ;
        %close all
 end






