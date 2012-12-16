function [ceq] = SolveImpCons(c1,R,x,s,Para)
% This function solves the implemetability conndition for c1 given x,R. We
% first express c2,l1,l2 in c1 and then use the result that
% x=xprime

%   Detailed explanation goes here
x=x;
n1=Para.n1;
n2=Para.n2;
g=Para.g(s);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma=Para.sigma;
c2=R^(-1/sigma)*c1;
% Solve for l1 , l2 using the resource constraint and wage equation
TotalResources=(c1*n1+c2*n2+g);
DenL2=theta_2*R*n1+theta_2*n2;
l2=(TotalResources-theta_1*n1+ theta_2*n1*R)/(DenL2);
if theta_2==0
    l1=TotalResources/(n1*theta_1);
l2=0;
else
l1= 1-(1-l2)*theta_2/theta_1*R;
end
% Nonlinear equality constraints - Imp Cons
ceq=(c2-c1)*(psi/(c2^sigma)) +((l1/(1-l1))*R-l2/(1-l2))*(1-psi)+x-x/beta;
end

