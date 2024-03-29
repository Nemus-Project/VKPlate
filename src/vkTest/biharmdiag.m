function [d0,d1,dc,dM1,dM] = biharmdiag(BCs,h,D,nu,dia) 
%BIHARMDIAG Summary of this function goes here
%   Detailed explanation goes here

pack_BCs = num2cell(BCs);
[K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};

[D00u00,~,~,~,~,~] = D00_coeffs(K0y,R0y,Kx0,Rx0,h,D,nu) ;
[D01u01,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,Rx0,h,D,nu) ;
[D02u02,~,~,~,~,~,~,~,~] = D02_coeffs(K0y,R0y,h,D,nu);
[D0Nu0N,~,~,~,~,~] = D00_coeffs(K0y,R0y,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,~,~,~,~,~,~,~] = D01_coeffs(K0y,R0y,RxL,h,D,nu) ;

[D10u10,~,~,~,~,~,~,~]           = D10_coeffs(R0y,Kx0,Rx0,h,D,nu) ;
[D11u11,~,~,~,~,~,~,~,~,~,~]     = D11_coeffs(R0y,Rx0,h,D,nu) ;
[D12u12,~,~,~,~,~,~,~,~,~,~,~]   = D12_coeffs(R0y,h,D,nu) ;
[D1Nu1N,~,~,~,~,~,~,~]           = D10_coeffs(R0y,KxL,RxL,h,D,nu) ;
[D1Nm1u1Nm1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(R0y,RxL,h,D,nu) ;

[D20u20,~,~,~,~,~,~,~,~]           = D20_coeffs(Kx0,Rx0,h,D,nu) ;
[D21u21,~,~,~,~,~,~,~,~,~,~,~]     = D21_coeffs(Rx0,h,D,nu) ;
[~,~,~,~,~,~,D22u22,~,~,~,~,~,~]   = D22_coeffs ;
[D2Nu2N,~,~,~,~,~,~,~,~]           = D20_coeffs(KxL,RxL,h,D,nu) ;
[D2Nm1u2Nm1,~,~,~,~,~,~,~,~,~,~,~] = D21_coeffs(RxL,h,D,nu) ;

[D10u10,~,~,~,~,~,~,~]    = D10_coeffs(RLy,Kx0,Rx0,h,D,nu) ;
[D11u11,~,~,~,~,~,~,~,~,~,~]     = D11_coeffs(RLy,Rx0,h,D,nu) ;
[D12u12,~,~,~,~,~,~,~,~,~,~,~]   = D12_coeffs(RLy,h,D,nu) ;
[D1Nm1u1Nm1,~,~,~,~,~,~,~,~,~,~] = D11_coeffs(RLy,RxL,h,D,nu) ;
[D1Nu1N,~,~,~,~,~,~,~]    = D10_coeffs(RLy,KxL,RxL,h,D,nu) ;

[D00u00,~,~,~,~,~] = D00_coeffs(KLy,RLy,Kx0,Rx0,h,D,nu) ;
[D01u01,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,Rx0,h,D,nu) ;
[D02u02,~,~,~,~,~,~,~,~] = D02_coeffs(KLy,RLy,h,D,nu) ;
[D0Nu0N,~,~,~,~,~] = D00_coeffs(KLy,RLy,KxL,RxL,h,D,nu) ;
[D0Nm1u0Nm1,~,~,~,~,~,~,~] = D01_coeffs(KLy,RLy,RxL,h,D,nu) ;

d0 = [D00u00; D01u01; D02u02; D0Nm1u0Nm1; D0Nu0N];
d1 = [D10u10; D11u11; D12u12; D1Nm1u1Nm1; D1Nu1N];
dc = [D20u20; D21u21; D22u22; D2Nm1u2Nm1; D2Nu2N];
dM1 = [D10u10; D11u11; D12u12; D1Nm1u1Nm1; D1Nu1N];
dM = [D00u00; D01u01; D02u02; D0Nm1u0Nm1; D0Nu0N];

D0 = [d0,d1,dc,dM1,dM];

end

