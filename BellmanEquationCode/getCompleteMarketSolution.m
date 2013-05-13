function [ PolicyRule,VNew,zLS,ifail ] = getCompleteMarketSolution(x,R,s_,Para,z0)
%getCompletMarketSoution ; This functioon solves the comlete market soltion
%associated with x, R, s_
    S = length(Para.P(s_,:));
    psi=Para.psi;
    sigma=Para.sigma;
    beta=Para.beta;
    P=Para.P;
    %z0 = [0.5*ones(1,4*S) zeros(1,3*S+1)];
    partialCompleteMarketResidual=@(n,z,user,iflag) completeMarketsResiduals(n,z,user,iflag,x,R,s_,Para);
    [zLS, ~,~,ifail]=c05qb(partialCompleteMarketResidual ,z0,'xtol',1e-10);
    c1 = zLS(1:S);
    c2 = zLS(S+1:2*S);
    l1 = zLS(2*S+1:3*S);
    l2 = zLS(3*S+1:4*S);
    xprime = x *ones(1,S);
    u = Para.alpha_1*uAlt(c1,l1,psi,sigma)+Para.alpha_2*uAlt(c2,l2,psi,sigma);
    A=eye(S)-(repmat(beta,S,1).*P);
    uu=(A\u')';
    VNew = sum(Para.P(s_,:).*uu);
    btildprime = Para.beta(1,:).*xprime./(psi*c2(1,:).^(-sigma));
    PolicyRule=[c1 c2 l1 l2 btildprime R.*ones(1,S) xprime];

end

