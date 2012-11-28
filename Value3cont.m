function [ minusVobj,minusGrad] = Value3cont(x)
global V Vcoef R u2btild Par s_
%BELOBJECTIVEUNCOND Computes the Bellman objective with 
%   Detailed explanation goes here
    psi = Par.psi;
    sigma = Par.sigma;
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
    [c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [Rprime,gradRprime] = computeR(c1,c2,gradc1,gradc2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,u2btild);
    
    
    
    
 
 
    %compute objective
    Vprime(:,1) = funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)])*ones(3,1);
    Vprime(:,2) = funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)])*ones(3,1);
    V_x(:,1)=funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)],[1,0])*ones(3,1);
    V_x(:,2)=funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)],[1,0])*ones(3,1);
    V_R(:,1)=funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)],[0,1])*ones(3,1);
    V_R(:,2)=funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)],[0,1])*ones(3,1);
    
    Vrhs = alpha(1)*uAlt(c1,l1,psi,sigma)+alpha(2)*uAlt(c2,l2,psi,sigma) + beta*Vprime;
    

    
    
gradV=alpha(1).*psi.* c1.^(-sigma).*gradc1...
        +alpha(2).*psi.* c2.^(-sigma).*gradc2...
        -alpha(1).*(1-psi)./(1-l1).*gradl1...
        -alpha(2).*(1-psi)./(1-l2).*gradl2...
        +beta*(V_x.*gradxprime+V_R.*gradRprime);
    
    
    minusGrad =-gradV*P(s_,:)';
    minusVobj = -Vrhs(1,:)*P(s_,:)';
     if max([l1 l2]) >1
                grad=abs(x)+100;

    end
    if ~isreal(minusGrad)
    
    minusGrad=-abs(minusGrad)-100;
    end
    else
        minusGrad=-abs(x)-100;
   
    
    end
     

end