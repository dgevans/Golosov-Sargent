function [res]=FOCFB(x,s,Para)
% FOC for the FB
psi=Para.psi;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1;
theta_2=Para.theta_2;
n1=Para.n1;
n2=Para.n2;
g=Para.g(s);
if min(x)>0
c1=x(1);
c2=x(2);
l1=x(3);
l2=x(4);
res(1)=alpha_1*c1^(-1)-alpha_2*c2^(-1);
res(2)=alpha_1/(theta_1*(1-l1))-alpha_2/(theta_2*(1-l2));
res(3)=(1-psi)/(theta_2*(1-l2))-psi/c2;
res(4)=c1*n1+c2*n2+g-theta_1*l1*n1-theta_2*l2*n2;
else
    res=abs(x)*10+100;
end
