function [resINQ,resEQ]=getResInitialLaborKnitro(x,s_,n,Para)
resINQ=[];
g=Para.g;
psi=Para.psi;
pi=Para.pi;
beta=Para.beta;
sSize=Para.sSize;
c=n-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
resEQ=zeros(sSize,1);
%if min(n-g)>0 && min(1-n)>0
for s=1:sSize
  resEQ(s)=  x*(1-R(s))+psi-(1-psi)*((n(s))/(1-n(s)));
end
%else
 %   resINQ=(abs(n-g)+abs(1-n)).*10;
end
