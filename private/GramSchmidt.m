function [PhiGS] = GramSchmidt(Phi,Nmodes)

PhiGS=zeros(length(Phi),Nmodes)
PhiGS(:,1)=Phi(:,1);
add=zeros(length(Phi),1);
for i= 2:Nmodes

   %preproj= preproj+;
    for j =1:i-1
        add=add+(dot(PhiGS(:,j),Phi(:,i))/(dot(PhiGS(:,j),PhiGS(:,j))))*PhiGS(:,j);
    end

   PhiGS(:,i)=Phi(:,i)-add;

end

% for i= 2:Nmodes
%     R = norm(PhiGS(:,i));
%   PhiGS(:,i)=PhiGS(:,i)/R
% end
