function [ v] = getValue0(x, btild,s_,Para,c,V)
% This function gets the time 0 value at the guess given by c1 and c2 for a
% given s0=s_


% RETRIVE THE PARAMETERS FROM THE STRUCT
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
g=Para.g(s_);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma=Para.sigma;
c1 = x(1);
c2 = x(2);
% CHECK IF CONSUMPTION IS NON NEGATIVE
if min(x)>0
% USE THE TIME0 CONSTRAINTS TO FIGURE OUT L1,L2,X0,R0    
    R=(c1/c2)^(sigma);
TotalResources=(c1*n1+c2*n2+g);
DenL2=theta_2(s_)*R*n1+theta_2(s_)*n2;
l2=(TotalResources-theta_1(s_)*n1+ theta_2(s_)*n1*R)./(DenL2);
l1= 1-(1-l2).*theta_2(s_)./theta_1(s_)*R;
xprime=-(c2-c1)*(psi*c2^(-sigma))-((l1/(1-l1)).*R-l2./(1-l2))*(1-psi)+btild*psi*c2^(-sigma);
Rprime=c2^(-sigma)/c1^(-sigma);
% COMPUTE THE VALUE AT TIME 0 USING THE POLICIES AND CONTINUATION VALUES AT
% T=1 THRU C,V
v=alpha_1*uAlt(c1,l1,psi,sigma)+alpha_2*uAlt(c2,l2,psi,sigma)+beta(s_)*funeval(c(s_,:)',V(s_),[xprime Rprime]);
v=-v;
else
    v=(max(abs(x)))*100+10;
    
end

