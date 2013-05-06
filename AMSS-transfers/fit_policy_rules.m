function [N coeffN XPrime coeffXPrime] =fit_policy_rules(s_,coeff,V,domain, nstore,Para,gridfactor,flagtr)

ApproxMethod='spli';
OrderOfApprx=25;
orderspli=3;

xGridSize=length(Para.xGrid)*gridfactor;
xGrid=linspace(Para.xMin,Para.xMax,xGridSize);

for xind=1:xGridSize
    x=xGrid(xind);
[nguess,~]=GetInitialApproxPolicy(x,domain,squeeze(nstore(:,s_,:)));  
[n,xprime,c,~,exitflag] =solveInnerOpt(x,s_,coeff,V,Para,nguess);
if flagtr==1

    [n,xprime,c,~,~] =solveInnerOptUsingConsOptFminCon(x,s_,coeff,V,Para,[n,xprime]); 
end
n_fit(xind,:)=n;
xprime_fit(xind,:)=xprime;

end


for s=1:sSize
N(s) = fundefn(ApproxMethod,OrderOfApprx ,Para.xMin,Para.xMax,orderspli);
XPrime(s)=fundefn`(ApproxMethod,OrderOfApprx ,Para.xMin,Para.xMax,orderspli);
coeffN(s,:)=funfitxy(N(s),xGrid,n_fit(:,1,s));
coeffXPrime(s,:)=funfitxy(XPrime(s),xGrid,xprime_fit(:,1,s));
end

end
