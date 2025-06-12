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
        % if nref==2
        %     if nQ==2
        % figure
        % pcolor(X,Y,(Spreftemp))
        % 
        %     end
        % end

        Qinterp= reshape(Spreftemp,[(Ny+1)*(Nx+1),1]) ;

        %trapzIntcalc(Qinterp.*Qinterp,h,Nx,Ny)

        Phiref(:,nQ) = Qinterp;
        Phicur(:,nQ) = Qtemp;
       
      
        MacMat(nref,nQ) = (abs(Phiref(:,nQ)'*Phicur(:,nQ)))^2/((Phiref(:,nQ)'*Phiref(:,nQ))*(Phicur(:,nQ)'*Phicur(:,nQ)));

      
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
        check=row(loo);
    end
end

% Phireg=ones((Ny+1)*(Nx+1),2);
% for nfix = 1 : Nmodes
%     if nfix==2
%         % A=((Phicur(:,nfix)'+ Phicur(:,nfix+1)')*Phiref(:,nfix))/(2*Phiref(:,nfix)'*Phiref(:,nfix));
%         % B=((Phicur(:,nfix)'- Phicur(:,nfix+1)')*Phiref(:,nfix+1))/(2*Phiref(:,nfix+1)'*Phiref(:,nfix+1));
% 
%         %A=Phicur(:,nfix)'*Phiref(:,nfix)
%         %B=Phicur(:,nfix+1)'*Phicur(:,nfix)
%         % C=Phicur(:,nfix)'*Phiref(:,nfix+1)
%         % D=Phicur(:,nfix+1)'*Phiref(:,nfix+1)
% 
%         %regMAT=inv([A,B;A,-B]);
%         %trapzIntcalc((Phicur(:,nfix)+ Phicur(:,nfix+1)).*(Phicur(:,nfix)+ Phicur(:,nfix+1)),h,Nx,Ny)
%        %((Phicur(:,nfix)'+ Phicur(:,nfix+1)')*Phiref(:,nfix))
% 
%         %Phireg=regMAT*Phicur(:,[nfix,nfix+1])';
%     end 
% end
% % return
% QMAC(:,[2,3])=Phireg';
% %pause(2)





