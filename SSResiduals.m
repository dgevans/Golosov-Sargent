function [ res ] = SSResiduals( X,Para )
%SSRESIDUALS Summary of this function goes here
%   Detailed explanation goes here
psi = Para.psi;
beta = Para.beta;
alpha_1 = Para.alpha_1;
alpha_2 = Para.alpha_2;
theta_1 = Para.theta_1;
theta_2 = Para.theta_2;
n1=Para.n1;
n2=Para.n2;
g = Para.g;
sigma = Para.sigma;
P = Para.P;
P = P(1,:);
Palt = fliplr(P);

c1 = X(1:2);
c2 = X(3:4);
l1 = X(5:6);
l2 = X(7:8);
R = X(9);
x = X(10);
mu = X(11);
lambda = X(12);
phi = X(13:14);
xi = X(15:16);
rho = X(17:18);

uc1 = psi*c1.^(-sigma);
uc2 = psi*c2.^(-sigma);
ul1 = -(1-psi)./(1-l1);
ul2 = -(1-psi)./(1-l2);
uc2Alt = fliplr(uc2);

Euc2 = dot(P,uc2);
Euc1 = dot(P,uc1);

res(1:2) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - x +R*l1.*ul1-l2.*ul2;

res(3:4) = theta_2*R*(1-l2)-theta_1*(1-l1);

res(5:6) = n1*theta_1*l1+n2*theta_2*l2-n1*c1-n2*c2-g;

res(7:8) = uc1*R-uc2;

res(9:10) = alpha_1*P.*uc1+mu*P.*uc2-n1*xi - sigma*R*rho.*uc1./c1;

res(11:12) = alpha_2*P.*uc2+mu*P.*uc2.*( (sigma*x*(P.*uc2-Euc2))./(beta*c2*Euc2^2)...
    +sigma-1-sigma*c1./c2 ) + (sigma*x*mu*Palt.*uc2Alt.*P.*uc2)./(beta*c2*Euc2^2)...
    -n2*xi+sigma*rho.*uc2./c2;

res(13:14) = alpha_1*P.*ul1-mu*(1-psi)*R*P./( (1-l1).^2 ) + (phi+n1*xi)*theta_1;

res(15:16) = alpha_2*P.*ul2+mu*(1-psi)*P./( (1-l2).^2 )+(n2*xi-R*phi)*theta_2;

res(17:18) = -beta*lambda*P*Euc1+mu*P.*ul1.*l1+lambda*P.*uc1+theta_2*phi.*(1-l2)...
    +rho.*uc1;




% res(1:2) =x*uc2/(beta*Euc2)-( psi*uc2.*(c2-c1)+x+(1-psi)*( R*l1./(1-l1) - l2./(1-l2) ) );
% 
% res(3:4) = theta_1*(1-l1)-theta_2*(1-l2)*R;
% 
% res(5:6) = theta_1*l1+theta_2*l2-c1-c2-g;
% 
% res(7:8) = psi*alpha_1*P.*uc1-psi*mu*P.*uc2-xi;
% 
% res(9:10) = alpha_2*psi*P.*uc2-mu*P.*( psi*(sigma*c1-(sigma-1)*c2)./(c2.^(sigma+1)) ...
%     -sigma*x*(Euc2-P.*uc2)./(c2.^(sigma+1)*Euc2^2) )+sigma*x*mu*Palt.*P.*uc2Alt./(beta*c2.^(sigma+1).*Euc2^2) ;
% 
% res(11:12) = theta_1*(xi+phi)-alpha_1*(1-psi)*P./(1-l1) -(1-psi)*R*mu*P./( (1-l1).^2 );
% 
% res(13:14) = theta_2*(xi-R*phi) - alpha_2*(1-psi)*P./(1-l2) + mu.*P./( (1-l2).^2 );
% 
% res(15:16) = beta*P*lambda*Euc1 - mu*(1-psi).*P.*l1./(1-l1)...
%     -lambda*P.*uc1+phi.*(1-l2);


end

