function [ v] = getValue0(x, btild,s_,Para,c,V)
%Solves the value for a given value of c1, c2 and btild
%   Detailed explanation goes here
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
g=Para.g(s_);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
c1 = x(1);
c2 = x(2);
if min(x)>0
    R=c1/c2;
TotalResources=(c1*n1+c2*n2+g);
FF=R*theta_2/theta_1;
DenL2=n1*theta_1*FF+theta_2*n2;


l2=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
l1= 1-FF*(1-l2);
if theta_2==0
    l2=0;
    l1=TotalResources/theta_1;
end
BracketTerm=l2/(1-l2)-(l1/(1-l1))*R;

% this is the IMPLEMENTABILITY after dividing be c2
u2btildprime=(((1-psi)/(psi))*BracketTerm+btild/(beta*psi)+R-1)*psi; % <CHECK THIS>
btildprime=u2btildprime/(c2^-1*psi) ;
Rprime=c2^(-1)/c1^(-1);
v=alpha_1*uBGP(c1,l1,psi)+alpha_2*uBGP(c2,l2,psi)+beta*funeval(c(1,:)',V(1),[u2btildprime Rprime]);
v=-v;
else
    v=(max(abs(x)))*100+10;
    
end

