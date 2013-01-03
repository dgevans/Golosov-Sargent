function [ res ] = StationaryPolicyRes(  X,R,Para )
%FOCRESIDUALS Summary of this function goes here
%   Detailed explanation goes here
psi = Para.psi;
beta = Para.beta;
theta_1 = Para.theta_1;
theta_2 = Para.theta_2;
n1=Para.n1;
n2=Para.n2;
g = Para.g;
sigma = Para.sigma;
P = Para.P;
P = P(1,:);

c1 = X(1:3);
c2 = X(4:6);
l1 = X(7:9);
l2 = X(10:12);

uc1 = psi*c1.^(-sigma);
uc2 = psi*c2.^(-sigma);
ul1 = -(1-psi)./(1-l1);
ul2 = -(1-psi)./(1-l2);

Euc2 = dot(P,uc2);



res(1:3) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - x +R.*l1.*ul1-l2.*ul2;

res(4) = dot(P,uc2)./dot(P,uc2)-R;

res(4:5) = theta_2*R.*(1-l2)-theta_1*(1-l1);

res(6:7) = n1*theta_1*l1+n2*theta_2*l2-n1*c1-n2*c2-g;

res(8:9) = uc1.*R-uc2;

end

