function [ grad, user,iflag] = BelObjectiveUncondGradNAGBGP(n,x,user,iflag)
global V Vcoef R u2btild Par s_
%BELOBJECTIVEUNCOND Computes the Bellman objective with 
%   Detailed explanation goes here
    psi = Par.psi;
    beta =  Par.beta;
    P = Par.P;
    theta_1 = Par.theta(:,1);
    theta_2 = Par.theta(:,2);
    g = Par.g;
    alpha = Par.alpha;
    n1 = Par.n1;
    n2 = Par.n2;
    sigma = 1;
    
    frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
    
    if (min(x)>0 && frac>0)
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
    %compute components from unconstrained guess
    [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);

    [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);

    [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,psi,beta,P);

 
    %compute objective
    grad1 = zeros(3,1);
    X = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
    Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi)...
            +beta*funeval(Vcoef{1},V(1),X));
    dV = funeval(Vcoef{1},V(1),X,eye(2));
      % Direct gardients with c_1_1,c_1_2,c_2_!
    grad1(1) = P(s_,1)*(alpha(1)*psi*c1_1^(-1)+beta*c2_1^(-1)*dV(2)); %<ok - Anmol>
    grad1(2) = 0; %<ok - Anmol>
    grad1(3) = P(s_,1)*(alpha(2)*psi*c2_1^(-1)-beta*c2_1^(-2)*(psi*btildprime(1)*dV(1)+c1_1*dV(2))); %<ok - Anmol>
    
    grad1 = grad1+P(s_,1)*( psi*grad_btildprime(:,1)*c2_1^(-1)*beta*dV(1) ...
        - alpha(1)*(1-psi)*l1grad(:,1)/(1-l1(1))...
        -alpha(2)*(1-psi)*l2grad(:,1)/(1-l2(1))); %<ok - Anmol>
    
    grad2 = zeros(3,1);
    X = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
    Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi)...
            +beta*funeval(Vcoef{2},V(2),X));
    dV = funeval(Vcoef{2},V(2),X,eye(2));
    grad2(1) = 0; %<ok - Anmol>
    grad2(2) = P(s_,2)*(alpha(1)*psi*c1_2^(-1)+beta*c2_2^(-1)*dV(2)); %<ok - Anmol>
    grad2(3) = 0;
    
    d_c2_2 = P(s_,2)*(alpha(2)*psi*c2_2^(-1)-beta*c2_2^(-2)*(dV(2)*c1_2 + psi*btildprime(2)*dV(1)));%<ok - Anmol>
    
    grad2 = grad2 + d_c2_2*grad_c2_2;
    grad2 = grad2 + P(s_,2)*( grad_btildprime(:,2)*psi*c2_2^(-1)*beta*dV(1)...
        -alpha(1)*(1-psi)*l1grad(:,2)/(1-l1(2))...
        -alpha(2)*(1-psi)*l2grad(:,2)/(1-l2(2)));
    
    grad = grad1+grad2;
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
