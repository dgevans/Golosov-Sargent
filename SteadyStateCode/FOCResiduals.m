function [ res ] = FOCResiduals(  X,R,x,VR,Vx,Para )
%FOCRESIDUALS Summary of this function goes here
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

c1 = X(1:2);
c2 = X(3:4);
l1 = X(5:6);
l2 = X(7:8);
Rprime = X(9:10);
xprime = X(11:12);
mu = X(13:14);
lambda = X(15);
phi = X(16:17);
xi = X(18:19);
rho = X(20:21);

uc1 = psi*c1.^(-sigma);
uc2 = psi*c2.^(-sigma);
ul1 = -(1-psi)./(1-l1);
ul2 = -(1-psi)./(1-l2);
uc2Alt = fliplr(uc2);
muAlt = fliplr(mu);

Euc2 = dot(P,uc2);



res(1:2) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - xprime +Rprime.*l1.*ul1-l2.*ul2;

res(3) = dot(P,uc1.*(Rprime-R));

res(4:5) = theta_2*Rprime.*(1-l2)-theta_1*(1-l1);

res(6:7) = n1*theta_1*l1+n2*theta_2*l2-n1*c1-n2*c2-g;

res(8:9) = uc1.*Rprime-uc2;

res(10:11) = alpha_1*P.*uc1+mu.*uc2-n1*xi - sigma*rho.*Rprime.*uc1./c1...
    -lambda*sigma*P.*uc1.*(Rprime-R)./c1;

res(12:13) = alpha_2*P.*uc2+mu.*uc2.*( (sigma*x*(P.*uc2-Euc2))./(beta*c2*Euc2^2)...
    +sigma-1-sigma*c1./c2 ) + (sigma*x*muAlt.*uc2Alt.*P.*uc2)./(beta*c2*Euc2^2)...
    -n2*xi+sigma*rho.*uc2./c2;

res(14:15) = alpha_1*P.*ul1-(1-psi)*Rprime.*mu./( (1-l1).^2 ) + (phi+n1*xi)*theta_1;

res(16:17) = alpha_2*P.*ul2+(1-psi)*mu./( (1-l2).^2 )+(n2*xi-Rprime.*phi)*theta_2;

res(18:19) = beta*P.*Vx - mu;

res(20:21) = beta*P.*VR+mu.*ul1.*l1+lambda*P.*uc1+theta_2*phi.*(1-l2)...
    +rho.*uc1;

end

