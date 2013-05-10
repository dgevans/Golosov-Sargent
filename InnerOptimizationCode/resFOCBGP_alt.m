function [ res, user,iflag] = resFOCBGP_alt(n,z,user,iflag)
global V Vcoef R x Par s_ flagCons upperFlags lowerFlags
% Computes the gradient of the bellamn equation objective under the
%constraints that xLL <= xprime <= xUL.  Will follow most
%of the methodology as BelObjectiveUncondGradNAGBGP but with a few
%differences.  One of the choice variables will be xprime, with the
%constraint that xprime computed from c_1(1) c_1(2) c_2(1) is equal to
%xprime.

    %get parameters from par
    psi = Par.psi;
    beta =  Par.beta;
    P = Par.P;
    theta_1 = Par.theta_1;
    theta_2 = Par.theta_2;
    g = Par.g;
    alpha = Par.alpha;
    n1 = Par.n1;
    n2 = Par.n2;
    xLL=Par.xLL;
    xUL=Par.xUL;
    sigma = Par.sigma;
    S = length(P(1,:));
    z = z(:)';
      if numel(beta)==1
    beta=ones(2*S-1,S)*beta;
    else
    beta=repmat(beta,2*S-1,1);
    end
    
    %get c_1 and c_2 from z
    c1 = z(1:S);
    c2_ = z(S+1:2*S-1);
    P_ = P(s_,:); P_(S) = [];
    frac = (R*dot(P(s_,:),c1.^(-sigma)) - dot(P_,c2_.^(-sigma)))/P(s_,S);
    %it is necessary for frac to be positive in order for the constraints
    %to be satisfied.
    if (min(z(1:2*S-1))>0 && frac>0)
        innerFlags = 1-upperFlags-lowerFlags;
    %MuL and MuH are the Lagrange multipliers for the upper and lower
    %constraints
     MuL = lowerFlags.*z(2*S:3*S-1);
     MuH = upperFlags.*z(2*S:3*S-1);
     xprime = innerFlags.*z(2*S:3*S-1)+lowerFlags.*xLL+upperFlags.*xUL;
    
    %lambda_I are the lagrange multipliers that insure that  xprime =
    %xprime computed below.  This is a bit hockey but it works well
    lambda_I = z(3*S:4*S-1);
    res = zeros(4*S-1,1);
    
    %Compute choice variables and gradients using unconstrained code
    [c1,c2,gradc1,gradc2] = computeC2_2(c1,c2_,R,s_,P,sigma);
    [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
                                        
    [ xprimeMat,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,x);
    
    %compute grad of the objective using xprime
    V_x = zeros(1,S);
    V_R = zeros(2*S-1,S);
    %need for loop because there may be different coefficients
    for s = 1:S
        V_x(1,s)=funeval(Vcoef{s},V(s),[xprime(1,s) Rprime(1,s)],[1,0]);
        V_R(:,s)=funeval(Vcoef{s},V(s),[xprime(1,s) Rprime(1,s)],[0,1])*ones(2*S-1,1);
    end


    %put lambda into standard 2*S-1xS format
    lambda = kron(ones(2*S-1,1),lambda_I);
    
    %first compute gradient of the objective with respecto the first 3
    %variables c_1(1), c_1(2) and c_2(1).  Note these don't effect V_x
    %(since it depends only on xprime) so that term is not included.
    %Gradients are computed for each state
    
    gradV=alpha(1).*psi.* c1.^(-sigma).*gradc1...
        +alpha(2).*psi.* c2.^(-sigma).*gradc2...
        -alpha(1).*(1-psi)./(1-l1).*gradl1...
        -alpha(2).*(1-psi)./(1-l2).*gradl2...
        +beta.*(V_R.*gradRprime)...
        -lambda.*gradxprime;
    
    %Combine states using transition matrix
    grad =gradV*P(s_,:)';
     
    res(1:2*S-1)=grad;
    
    % Next we have the two first order conditions with respect to
    % xprime.
    %res(4) = P(s_,1)*lambda_I(1)+P(s_,1)*beta*V_x(1,1)+MuL(1)-MuH(1); 
    %res(5) = P(s_,2)*lambda_I(2)+P(s_,2)*beta*V_x(1,2)+MuL(2)-MuH(2); 
    res(2*S:3*S-1) = (P(s_,:).*(lambda_I+beta(1,:).*V_x)+MuL-MuH)';
   
    % FOC with respect to labmda_I imposing that xprime = xprim2
    %res(6) = xprime(1)-xprimeMat(1,1);
    %res(7) = xprime(2)-xprimeMat(1,2);
    res(3*S:4*S-1) = (xprime-xprimeMat(1,:))';
     
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
