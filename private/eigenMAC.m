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
       
        Phiref(:,nref)=Phi1;
        Phicur(:,nQ)=Phi2;

      
        MacMat(nref,nQ) = (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));

      
    end
end

MacMat(MacMat<0.1)=0;
for n = 1 : Nmodes
    MacMat(n,n)=0;
end

[row,col]=find(MacMat);
lmc=length(col);
%lis=[row,col];
if lmc~=0
    check=row(1);
end
for loo= 1:lmc-1
    if check~=col(loo)

        QMAC(:,[col(loo),row(loo)])=QMAC(:,[row(loo),col(loo)]);
        OmMAC([col(loo),row(loo)])=OmMAC([row(loo),col(loo)]);
        check=row(loo);
    end
end

%%
Omcheck=[Om(2:end);0];
degencheck=1./abs(Omcheck-Om);
degenum=find(abs(degencheck(:))>1);
%degenum2=degenum+1;

Phireg=ones((Ny+1)*(Nx+1),2);
for nfix = 1 : Nmodes
    if ismember(nfix, degenum)
        A=((Phicur(:,nfix)'+ Phicur(:,nfix+1)')*Phiref(:,nfix))/(2*Phiref(:,nfix)'*Phiref(:,nfix));
        B=((Phicur(:,nfix)'- Phicur(:,nfix+1)')*Phiref(:,nfix+1))/(2*Phiref(:,nfix+1)'*Phiref(:,nfix+1));

        A=Phicur(:,nfix)'*Phiref(:,nfix);
        B=Phicur(:,nfix+1)'*Phiref(:,nfix);
        C=Phicur(:,nfix)'*Phiref(:,nfix+1);
        D=Phicur(:,nfix+1)'*Phiref(:,nfix+1);
 
         regMAT=[A,B;C,D];
         %trapzIntcalc((Phi2(:,nfix)+ Phi2(:,nfix+1)).*(Phi2(:,nfix)+ Phi2(:,nfix+1)),h,Nx,Ny)
        
         Phireg=regMAT*Phicur(:,[nfix,nfix+1])';
         QMAC(:,[nfix,nfix+1])=Phireg';
     end 
end

%modelist=degenum;

 %QMAC(:,[2,3])=Phireg';
 %QMAC(:,[2,3])=Phiref(:,[2,3]);



