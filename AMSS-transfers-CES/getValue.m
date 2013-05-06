function [val]=getValue(yy,x,s_,coeff,V,Para)

g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
n=yy(1:sSize);
xprime=yy(sSize+1:end);
c=n-g;

val=(psi*log(c(1))+(1-psi)*log(1-n(1))+beta*funeval(coeff(1,:)',V(1),xprime(1)))*pi(s_,1);
val=val+(psi*log(c(2))+(1-psi)*log(1-n(2))+beta*funeval(coeff(2,:)',V(2),xprime(2)))*pi(s_,2);
val=-val;


