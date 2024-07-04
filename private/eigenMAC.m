function [QMAC,OmMAC] = eigenMAC(Qref,Nxref,Nyref,Q,Nx,Ny,h,Nmodes,Lx,Ly,Om)

xref=linspace(0,Lx,Nxref+1);
yref=linspace(0,Ly,Nyref+1);
[XREF,YREF]=meshgrid(xref,yref);
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
[X,Y]=meshgrid(x,y);
h=1;
OmMAC=Om;
QMAC=Q;
MacMat = zeros(Nmodes,Nmodes);
for nQ = 1 : Nmodes
    for nref = 1 : Nmodes
        Qtemp   = Q(:,nQ) ;

        %Sptemp = reshape(Qtemp,[(Ny+1),(Nx+1)]) ;

        Qreftemp= Qref(:,nref) ;

        Spreftemp=reshape(Qreftemp,[(Nyref+1),(Nxref+1)]) ;

        Spreftemp= interp2(XREF,YREF,Spreftemp,X,Y);

        Qinterp= reshape(Spreftemp,[(Ny+1)*(Nx+1),1]) ;

        Phi1 = Qinterp;
        Phi2 = Qtemp;
       
      
        MacMat(nref,nQ) = (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));

      
    end
end

MacMat(MacMat<0.1)=0;
for n = 1 : Nmodes
    MacMat(n,n)=0;
end

[row,col]=find(MacMat)
lmc=length(col);
%lis=[row,col];
if lmc~=0
    check=row(1);
end
for loo= 1:lmc-1
    if check~=col(loo)

        QMAC(:,[col(loo),row(loo)])=QMAC(:,[row(loo),col(loo)]);
        OmMAC([col(loo),row(loo)])=OmMAC([row(loo),col(loo)]);
        check=row(loo)
    end
end

%pause(2)





