function [ceq] = ImpCons(x,R,u2btild,Para)
% Nonlinear inequality constraints
c = [];
%Solves the value for a given value of c1, c2 and btild
%   Detailed explanation goes here
c1 = x(1);
c2 = x(2);

n1=Para.n1;
n2=Para.n2;
g=Para.g(1);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
sigma=Para.sigma;
gamma=Para.gamma;
beta=Para.beta;

% Solve for l1 , l2 using the resource constraint and wage equation
TotalResources=(c1*n1+c2*n2+g);
DenL1=theta_1*n1+(c2/c1)^(-sigma/gamma)*(theta_2/theta_1)^(1/gamma)*n2*theta_2;
l1=TotalResources/DenL1;
l2=(c2/c1)^(-sigma/gamma)*(theta_2/theta_1)^(1/gamma)*l1;

% Nonlinear equality constraints - Imp Cons
% Check this - Anmol
ceq(1) = (n2*c2-n1*c1)-c2^(sigma)*n2*l2^(1+gamma)+n1*l1^(1+gamma)*c1^(sigma)-c2^(sigma)*u2btild*(1/beta-1);
ceq(2)=R-c2^(-sigma)/c1^(-sigma);

end

