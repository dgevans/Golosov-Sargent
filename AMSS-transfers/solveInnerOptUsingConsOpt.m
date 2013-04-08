function [n,xprime,c,VNew,exitflag] =solveInnerOptUsingConsOpt(x,s_,coeff,V,Para,z)
warning('on', 'NAG:warning')
g=Para.g;
psi=Para.psi;
beta=Para.beta;
xMax=Para.xMax;
xMin=Para.xMin;
options=Para.options;
pi=Para.pi;
sSize=Para.sSize;
nguess=[Para.nDet(x) Para.nDet(x)];
c=nguess-g;
R=1./(beta*(c).*sum(pi(s_,:).*(1./c)));
% Use Implementability to get xprime
xprimeguess=(nguess./(1-nguess))*(1-psi)+ x.*R-psi;

a = [];
bl = [Para.g(1); Para.g(2); Para.xMin; Para.xMin; -Inf;-Inf];
bu = [1; 1; Para.xMax; Para.xMax; 0;0];
istate = zeros(6, 1, 'int64');
cjac = zeros(2, 4);
clamda = zeros(6, 1);
r = zeros(4, 4);
%z = [nguess xprimeguess];
confun=@(mode, ncnln, n, ldcj, needc, z,ImpConsJac, nstate, user) ImplementabilityCons(mode, ncnln, n, ldcj, needc, z, x,s_,Para,ImpConsJac, nstate, user);
objfun=@(mode, n, z,objgrd, nstate, user) ValueFunction(mode, n, z,s_,coeff,V,Para,objgrd, nstate, user);
[cwsav,lwsav,iwsav,rwsav,ifail] = nag_opt_init('nag_opt_nlp1_solve');
[iter, istate, c, cjac, clamda, objf, objgrd, r, z] = ...
  nag_opt_nlp1_solve(a, bl, bu, confun, objfun, istate, cjac, clamda, r, z, lwsav, iwsav, rwsav);
n=z(1:2);
xprime=z(3:4);
c=n-g;
exitflag=ifail;
VNew=-objf;
end