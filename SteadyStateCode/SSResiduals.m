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
S = length(P);
c1 = X(1:S);
c2 = X(S+1:2*S);
l1 = X(2*S+1:3*S);
l2 = X(3*S+1:4*S);
R = X(4*S+1);
x = X(4*S+2);
mu = X(4*S+3);
lambda = X(4*S+4);
phi = X(4*S+5:5*S+4);
xi = X(5*S+5:6*S+4);
rho = X(6*S+5:7*S+4);

uc1 = psi*c1.^(-sigma);
uc2 = psi*c2.^(-sigma);
ul1 = -(1-psi)./(1-l1);
ul2 = -(1-psi)./(1-l2);

Euc2 = dot(P,uc2);
Euc1 = dot(P,uc1);

Euc2Grad = sum(sigma*x*mu*P.*uc2/(beta*Euc2^2));

res(1:S) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - x +R*l1.*ul1-l2.*ul2;

res(S+1:2*S) = theta_2*R*(1-l2)-theta_1*(1-l1);

res(2*S+1:3*S) = n1*theta_1*l1+n2*theta_2*l2-n1*c1-n2*c2-g;

res(3*S+1:4*S) = uc1*R-uc2;

res(4*S+1:5*S) = alpha_1*P.*uc1+mu*P.*uc2-n1*xi - sigma*R*rho.*uc1./c1;

res(5*S+1:6*S) = alpha_2*P.*uc2+mu*P.*uc2.*( (sigma*x*(-Euc2))./(beta*c2*Euc2^2)...
    +sigma-1-sigma*c1./c2 ) + Euc2Grad*P.*uc2./c2...
    -n2*xi+sigma*rho.*uc2./c2;

res(6*S+1:7*S) = alpha_1*P.*ul1-mu*(1-psi)*R*P./( (1-l1).^2 ) + (phi+n1*xi)*theta_1;

res(7*S+1:8*S) = alpha_2*P.*ul2+mu*(1-psi)*P./( (1-l2).^2 )+(n2*xi-R*phi)*theta_2;

res(8*S+1:9*S) = -beta*lambda*P*Euc1+mu*P.*ul1.*l1+lambda*P.*uc1+theta_2*phi.*(1-l2)...
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

