function [ v] = getValue0(x, btild,Para,c,V)
%Solves the value for a given value of c1, c2 and btild
%   Detailed explanation goes here
c1 = x(1);
c2 = x(2);
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
g=Para.g(1);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
sigma=Para.sigma;
gamma=Para.gamma;


TotalResoucrces=(c1*n1+c2*n2+g(1));
DenL1=theta_1*n1+(c2/c1)^(-sigma/gamma)*(theta_2/theta_1)^(1/gamma)*n2*theta_2;
l1=TotalResoucrces/DenL1;
l2=(c2/c1)^(-sigma/gamma)*(theta_2/theta_1)^(1/gamma)*l1;

btildprime=c2^(sigma)*btild/beta2-(c2-c1-l2^(1+gamma)*c2^(sigma)+l1^(1+gamma)*c1^(sigma));
u2btildprime=c2^(-sigma)*btildprime;
Rprime=c2^(-sigma)/c1^(-sigma);
v=alpha_1*u(c1,l1,sigma,gamma)+alpha_2*u(c2,l2,sigma,gamma)+beta*beta*funeval(c(1,:),V(1),[u2btildprime Rprime]);
v=-v;
end

