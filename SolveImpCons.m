function [ceq] = SolveImpCons(c1,R,u2btild,s,Para)
% This function solves the implemetability conndition for c1 given x,R. We
% first express c2,l1,l2 in c1 and then use the result that
% u2btild=u2btildprime

%   Detailed explanation goes here



n1=Para.n1;
n2=Para.n2;
g=Para.g(s);
theta_1=Para.theta_1(s);
theta_2=Para.theta_2(s);
psi=Para.psi;
beta=Para.beta;
c2=R^(-1)*c1;
% Solve for l1 , l2 using the resource constraint and wage equation
TotalResources=(c1*n1+c2*n2+g);
FF=R*theta_2/theta_1;
DenL2=n1*theta_1*FF+theta_2*n2;
l2=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
l1= 1-FF*(1-l2);
cUpperBound=(theta_1*n1*FF+theta_2*n2-theta_1*n1*(FF-1)-g)/(n1+n2/R);
% Nonlinear equality constraints - Imp Cons
%ceq = (c2-c1)-((1-psi)/psi)*((l2*c2)/(1-l2)-(l1*c1)/(1-l1))-(c2/psi)*u2btild*(1/beta-1);
BracketTerm=l2/(1-l2)-(l1/(1-l1))*R;
ceq = 1-R+u2btild/psi-((1-psi)/(psi))*BracketTerm-u2btild/(beta*psi);

end

