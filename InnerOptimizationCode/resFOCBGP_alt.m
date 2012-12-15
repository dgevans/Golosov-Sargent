function [ res, user,iflag] = resFOCBGP_alt(n,z,user,iflag)
global V Vcoef R u2btild Par s_ flagCons
%BELOBJECTIVEUNCOND Computes the gradient of the bellamn equation objective under the
%constraints that u2btildLL <= u2btildprime <= u2btildUL.  Will follow most
%of the methodology as BelObjectiveUncondGradNAGBGP but with a few
%differences.  One of the choice variables will be u2btildprime, with the
%constraint that xprime computed from c_1(1) c_1(2) c_2(1) is equal to
%u2btildprime.

    %get parameters from par
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
    sigma = Par.sigma;
    
    
    
    %Check if frac >0 to insure that x is a good guess.
    frac = (R*P(s_,1)*z(1)^(-sigma)+R*P(s_,2)*z(2)^(-sigma)-P(s_,1)*z(3)^(-sigma))...
        /( P(s_,2) );
    if (min(z(1:3))>0 && frac>0)
     %get c_1 and c_2 from z
     c1_1=z(1);
     c1_2=z(2);
     c2_1=z(3);
     
    %MuL and MuH are the Lagrange multipliers for the upper and lower
    %constraints
     MuL = zeros(2,1);
     MuH = zeros(2,1);
    switch flagCons
    case 'LL_'
       % lower limit binds for state 1 only.  So z(4) is the lagrange
       % multiplier for for lower constraint in state 1.  u2btildprime(2)
       % (also known as xprime) is allowed to move freely.  u2btildprime(1)
       % is set at the lower constraint
       MuL(1)=z(4);
       MuL(2)=0;
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=z(5);
       
    case '_LL'
       % lower limit binds for state 2 only. So z(5) is the lagrange
       % multiplier for for lower constraint in state 2.  u2btildprime(1)
       % (also known as xprime) is allowed to move freely.  u2btildprime(2)
       % is set at the lower constraint
       MuL(1)=0;
       MuL(2)=z(5);
       u2btildprime(1)=x(4);
       u2btildprime(2)=u2btildLL;
       
    case 'LLLL'
      % lower limit binds for both the states. z(4) and z(5) are the lower
      % limit Lagrange multipliers.  u2btildprime are set at the lower
      % limits
       MuL(1)=z(4);
       MuL(2)=z(5);
       u2btildprime(1)=u2btildLL;   
       u2btildprime(2)=u2btildLL;     
        
    case 'UL_'
     % upper limit binds for state 1 only.  So z(4) is the lagrange
       % multiplier for for upper constraint in state 1.  u2btildprime(2)
       % (also known as xprime) is allowed to move freely.  u2btildprime(1)
       % is set at the upper constraint

       MuH(1)=z(4);
       MuH(2)=0;
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=z(5);
       
        
    case '_UL'
         % upper limit binds for state 2 only  So z(5) is the lagrange
       % multiplier for for upper constraint in state 2.  u2btildprime(1)
       % (also known as xprime) is allowed to move freely.  u2btildprime(2)
       % is set at the upper constraint
       MuH(1)=0;
       MuH(2)=z(5);
       u2btildprime(1)=z(4);
       u2btildprime(2)=u2btildUL;
        
        
    case 'ULUL'
        
        
       % upper limit binds for both the states.  z(4) and z(5) are the
       % upper limit Lagrange multipliers.  u2btildprime are set at the
       % upper limits
       MuH(1)=z(4);
       MuH(2)=z(5);
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=u2btildUL;     
        
    otherwise
       MuL(1)=0;
       MuL(2)=0;
       u2btildprime(1)=z(4);
       u2btildprime(2)=z(5);     
        
    end
    
    %lambda_I are the lagrange multipliers that insure that  u2btildprime =
    %xprime computed below.  This is a bit hockey but it works well
    lambda_I(1) = x(6);
    lambda_I(2) = x(7);
    res = zeros(7,1);
    
    %Compute choice variables and gradients using unconstrained code
    [c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,u2btild);
    
    %compute grad of the objective using u2btildprime
    V_x(:,1)=funeval(Vcoef{1},V(1),[u2btildprime(1) Rprime(1,1)],[1,0])*ones(3,1);
    V_x(:,2)=funeval(Vcoef{2},V(2),[u2btildprime(2) Rprime(1,2)],[1,0])*ones(3,1);
    V_R(:,1)=funeval(Vcoef{1},V(1),[u2btildprime(1) Rprime(1,1)],[0,1])*ones(3,1);
    V_R(:,2)=funeval(Vcoef{2},V(2),[u2btildprime(2) Rprime(1,2)],[0,1])*ones(3,1);

    %put lambda into standard 3x2 format
    lambda = kron(ones(3,1),lambda_I);
    
    %first compute gradient of the objective with respecto the first 3
    %variables c_1(1), c_1(2) and c_2(1).  Note these don't effect V_x
    %(since it depends only on u2btildprime) so that term is not included.
    %Gradients are computed for each state
    gradV=alpha(1).*psi.* c1.^(-sigma).*gradc1...
        +alpha(2).*psi.* c2.^(-sigma).*gradc2...
        -alpha(1).*(1-psi)./(1-l1).*gradl1...
        -alpha(2).*(1-psi)./(1-l2).*gradl2...
        +beta*(V_R.*gradRprime)...
        -lambda.*gradxprime;
    
    %Combine states using transition matrix
    grad =gradV*P(s_,:)';
     
    res(1:3)=grad;
    
    % Next we have the two first order conditions with respect to
    % u2btildprime.
    res(4) = P(s_,1)*lambda_I(1)+P(s_,1)*beta*V_x(1,1)+MuL(1)-MuH(1); 
    res(5) = P(s_,2)*lambda_I(2)+P(s_,2)*beta*V_x(1,2)+MuL(2)-MuH(2); 
   
    % FOC with respect to labmda_I imposing that xprime = u2btildprim2
    res(6) = u2btildprime(1)-xprime(1,1);
    res(7) = u2btildprime(2)-xprime(1,2);

     
    %Various peices of code to return a large value if in a bad region.
     if max([l1(1,:) l2(1,:)]) >1
                res=abs(z)+100;

    end
     
    if ~isreal(grad)
    
    res=abs(z)+100;
    end
    else
        res=abs(z)+100;
    end

end
