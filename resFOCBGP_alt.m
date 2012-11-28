function [ res, user,iflag] = resFOCBGP_alt(n,x,user,iflag)
global V Vcoef R u2btild Par s_ flagCons
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
    u2btildLL=Par.u2btildLL;
    u2btildUL=Par.u2btildUL;
    sigma = 1;
    frac = (R*P(s_,1)*x(1)^(-1)+R*P(s_,2)*x(2)^(-1)-P(s_,1)*x(3)^(-1))...
        /( P(s_,2) );
    
    
    
    if (min(x(1:3))>0 && frac>0)
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
     MuL = zeros(2,1);
     MuH = zeros(2,1);
    switch flagCons
    case 'LL_'
       % lower limit binds for state 1 only
       MuL(1)=x(4);
       MuL(2)=0;
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=x(5);
       
    case '_LL'
       % lower limit binds for state 2 only
       MuL(1)=0;
       MuL(2)=x(5);
       u2btildprime(1)=x(4);
       u2btildprime(2)=u2btildLL;
       
    case 'LLLL'
      % lower limit binds for both the states
       MuL(1)=x(4);
       MuL(2)=x(5);
       u2btildprime(1)=u2btildLL;   
       u2btildprime(2)=u2btildLL;     
        
    case 'UL_'
     % upper limit binds for state 1 only

       MuH(1)=x(4);
       MuH(2)=0;
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=x(5);
       
        
    case '_UL'
         % upper limit binds for state 2 only
       MuH(1)=0;
       MuH(2)=x(5);
       u2btildprime(1)=x(4);
       u2btildprime(2)=u2btildUL;
        
        
    case 'ULUL'
        
        
       % upper limit binds for both the states
       MuH(1)=x(4);
       MuH(2)=x(5);
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=u2btildUL;     
        
    otherwise
       MuL(1)=0;
       MuL(2)=0;
       u2btildprime(1)=x(4);
       u2btildprime(2)=x(5);     
        
    end
    lambda_I(1) = x(6);
    lambda_I(2) = x(7);
    res = zeros(7,1);
    
    %compute components from unconstrained guess
   % [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,u2btild);
    
    %compute grad of the objective
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
     
    res(1:3)=grad;
    % FOC with respect to x'(1)
    res(4) = P(s_,1)*lambda_I(1)+P(s_,1)*beta*V_x(1,1)+MuL(1)-MuH(1); 
    res(5) = P(s_,2)*lambda_I(2)+P(s_,2)*beta*V_x(1,2)+MuL(2)-MuH(2); 
   
    % Definition of x'
    res(6) = u2btildprime(1)-xprime(1);
    res(7) = u2btildprime(2)-xprime(2);

           
     if max([l1(1,:) l2(1,:)]) >1
                res=abs(x)+100;

    end
     
    if ~isreal(grad)
    
    res=abs(x)+100;
    end
    else
        res=abs(x)+100;
    end

end
