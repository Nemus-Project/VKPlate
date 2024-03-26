function L = vkOperator(F1,F2,h,Nx,Ny)

Dx = DxBuild(Nx,Ny,h) ;
Dy = DyBuild(Nx,Ny,h) ;

L = (Dx*Dx*F1).*(Dy*Dy*F2) + (Dy*Dy*F1).*(Dx*Dx*F2) - 2*(Dx*Dy*F1).*(Dx*Dy*F2) ;

end
