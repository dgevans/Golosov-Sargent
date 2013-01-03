function [ res ] = StationaryPolicyRes( X,x,R,Para )
psi = Para.psi;
beta = Para.beta;
theta_1 = Para.theta_1;
theta_2 = Para.theta_2;
n1=Para.n1;
n2=Para.n2;
g = Para.g';
sigma = Para.sigma;
P = Para.P;
P = P(1,:);
n=length(g);
c1 = X(1:n);
c2 = X(n+1:2*n);
l1 = X(2*n+1:3*n);
l2 = X(3*n+1:4*n);

uc1 = psi*c1.^(-sigma);
uc2 = psi*c2.^(-sigma);
ul1 = -(1-psi)./(1-l1);
ul2 = -(1-psi)./(1-l2);

Euc2 = dot(P,uc2);



res(1:n) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - x +R.*l1.*ul1-l2.*ul2;
res(n+1:2*n) = theta_2*R.*(1-l2)-theta_1*(1-l1);
res(2*n+1:3*n) = n1*theta_1*l1+n2*theta_2*l2-n1*c1-n2*c2-g;
res(3*n+1) = dot(P,uc2)./dot(P,uc1)-R;

end

