function [res,user,iflag]=completeMarketsResiduals(n,z,user,iflag,x,R,s_,Para)
    %'''
    %Compute the residuals for the complete markets solution for a given state
    %'''
    %get some things from Para
    S = length(Para.P(s_,:));
    sigma=Para.sigma;
    beta = Para.beta;
    psi=Para.psi;
    P = Para.P;
    g = Para.g;
    n1 = Para.n1;
    n2 = Para.n2;
    alpha_1 =Para.alpha_1;
    alpha_2 = Para.alpha_2;
    theta_1 = Para.theta_1;
    theta_2 = Para.theta_2;
    c1 = z(1:S);
    c2 = z(S+1:2*S);
    l1 = z(2*S+1:3*S);
    l2 = z(3*S+1:4*S);
    phi = z(4*S+1:5*S);
    xi = z(5*S+1:6*S);
    rho = z(6*S+1:7*S);
    mu = z(7*S+1);
    checkc1=min(c1>0);
    checkc2=min(c2>0);
    checkl1=min(l1<1);
    checkl2=min(l2<1);
    if checkc1*checkc2*checkl1*checkl2==1
    uc1=psi*c1.^(-sigma);
    ucc1=psi*(-sigma).*c1.^(-sigma-1);
    ul1=-(1-psi)./(1-l1);
    ull1=-(1-psi)./((1-l1).^2);
     uc2=psi*c2.^(-sigma);
     ucc2=-psi*sigma*c2.^(-sigma-1);
     ul2=-(1-psi)./(1-l2);
      ull2=-(1-psi)./((1-l2).^2);
        
    res = zeros(7*S+1,1);
    
    res(1:S) = theta_2.*R.*ul1-ul2.*theta_1;
    res(S+1:2*S) = n1*l1.*theta_1+n2*l2.*theta_2 - n1*c1-n2*c2-g;
    res(2*S+1:3*S) = R*uc1-uc2;
    res(3*S+1:4*S) = alpha_1.*uc1+R.*mu.*( uc1+ucc1.*c1 ) -n1.*xi + ucc1.*R.*rho;
    res(4*S+1:5*S) = alpha_2.*uc2-mu.*( uc2+ucc2.*c2 ) - n2.*xi - ucc2.*rho;
    res(5*S+1:6*S) = alpha_1.*ul1 + R.*mu.*( ul1 + ull1.*l1 ) +theta_2.*R.*ull1.*phi./theta_1 + n1*theta_1.*xi;
    res(6*S+1:7*S) = alpha_2*ul2 - mu.*( ul2 + ull2.*l2 ) -ull2.*phi + n2*theta_2.*xi;
    I  = uc2.*c2+ul2.*l2-R*( uc1.*c1 + ul1.*l1);
    A=eye(S)-(repmat(beta,S,1).*P);
    xprime=(A\I')'; 
   res(7*S+1) = sum(P(s_,:).*xprime)-x;
    else
        res=abs([c1 c2 l1-1 l2-1 zeros(1,3*S) 1 ])+10;
    end
