function [ grad, user,iflag] = BelObjectiveUncondGradNAGBGP(n,x,user,iflag)
global V Vcoef R u2btild Par s_
%BELOBJECTIVEUNCOND Computes the Bellman objective with 
%   Detailed explanation goes here
    psi = Par.psi;
    beta =  Par.beta;
    P = Par.P;
    theta_1 = Par.theta(1);
    theta_2 = Par.theta(2);
    g = Par.g;
    alpha = Par.alpha;
    n1 = Par.n1;
    n2 = Par.n2;
    
    frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
    
    if (min(x)>0 && frac>0)
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
    %compute components from unconstrained guess
    [c1,c2,grad_c1,grad_c2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [Rprime,gradRprime] = computeR(c1,c2,grad_c1,grad_c2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_);
    
    
    
    
 
 
    %compute objective
    V_x(:,1)=funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)],[1,0])*ones(3,1);
    V_x(:,2)=funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)],[1,0])*ones(3,1);
    V_R(:,1)=funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)],[0,1])*ones(3,1);
    V_R(:,2)=funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)],[0,1])*ones(3,1);

    
    
gradV=alpha(1).*psi.* c1.^(-sigma).*gardc1...
        +alpha(2).*psi.* c2.^(-sigma).*gardc2...
        -alpha(1).*(1-psi)./(1-l1).*gradl1...
        -alpha(2).*(1-psi)./(1-l2).*gradl2...
        +beta*(V_x.*gradxprime+V_R.*gradRprime);
    
    
    grad =gradV.*P(s_,:)';
     if max([l1 l2]) >1
                grad=abs(x)+100;

    end
    if ~isreal(grad)
    
    grad=abs(grad)+100;
    end
    else
        grad=abs(x)+100;
   
    
    end
     

end