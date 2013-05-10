function [ grad, user,iflag] = BelObjectiveUncondGradNAGBGP(n,z,user,iflag)
global V Vcoef R x Par s_
%BELOBJECTIVEUNCOND Computes the gradient of the bellman equation objective with
%respect to c_1(1), c_1(2) and c_2(1).  Substitues out for the rest of
%the variables using their respective gradients.  z is the vector
%containing these three variables
%  Note that we pass global variables as well (because NAG did not allow us
%  to pass user defined information) V is a struct containing the
%  interpolation information.  Vcoef has coefficents for the value
%  function. R, x and s_ are the previous period states.  Par is a
%  struct containig all the relevent parameters.  We return the gradient
%  grad.  user and iflag our variables that nag requires but we don't use.

       %get parameters from the Par struct
    psi = Par.psi;
    sigma = Par.sigma;
    beta =  Par.beta;
    P = Par.P;
    theta_1 = Par.theta_1;
    theta_2 = Par.theta_2;
    g = Par.g;
    alpha = Par.alpha;
    n1 = Par.n1;
    n2 = Par.n2;
    
    %make sure z is a row vector
    z = z(:)';

    S = length(P(1,:));
    c1 = z(1:S);
    c2_ = z(S+1:2*S-1);
    P_ = P(s_,:); P_(S) = [];
    if numel(beta)==1
    beta=ones(2*S-1,S)*beta;
    else
    beta=repmat(beta,2*S-1,1);
    end
    frac = (R*dot(P(s_,:),c1.^(-sigma)) - dot(P_,c2_.^(-sigma)))/P(s_,S);
    %it is necessary for frac to be positive in order for the constraints
    %to be satisfied.
    if (min(z)>0 && frac>0)
    %compute components from unconstrained guess
    
   
    [c1,c2,gradc1,gradc2] = computeC2_2(c1,c2_,R,s_,P,sigma);
     %compute c_2(2) and return c1 and c2 as 3x2 matrix.  Here c1 will be of
    %the form
    %     c_1(1)   c_1(2)
    %     c_1(1)   c_1(2)
    %     c_1(1)   c_1(2)
    %Similarly for c2.  The 3x2 matrix format is useful for computing future
    %gradients.  gradc1 will a matrix containing the derivative of
    %c1 with respect to the various z's.  For instance the ith row and jth
    %column will be the derivative of c_1(j) with respect to x(i).  Thus
    %gradc1 will be
    %       1   0
    %       0   1      
    %       0   0
    %Similarly for gradc2
    
    
    [Rprime,gradRprime] = computeR(c1,c2,gradc1,gradc2,sigma);
    %Computes Rprime c_2(s)^(-sigma)/c_1(s)^(-sigma) in the 3x2 matrix form
    %described above as well as its gradient with respect to z
    
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2);
    %computes l_1 and l_2, the labor supply  of agent 1 and 2 in the
    %standard 3x2 format, along with their gradients with respect to z                                    
                                        
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,x);
    %Computes the choice of the state variable xprime tomorrow in the
    %standard 3x2 format as well as gradient with respect to z (note this
    %is unfortunated notation, xprime is refering to xprime
    %($u_{c,2}\tilde b'$), while z is the vector [c_1(1), c_1(2), c_2(1)].
    %Sorry for the confusion
    
    
    
    
 
 
    %compute objective derivative of the objectve with respect to the state
    %variables for both states (s = 1 or 2) of the world.
    
    V_x = zeros(2*S-1,S);
    V_R = zeros(2*S-1,S);
    %need for loop because there may be different coefficients
    for s = 1:S
        V_x(:,s)=funeval(Vcoef{s},V(s),[xprime(1,s) Rprime(1,s)],[1,0])*ones(2*S-1,1);
        V_R(:,s)=funeval(Vcoef{s},V(s),[xprime(1,s) Rprime(1,s)],[0,1])*ones(2*S-1,1);
    end
    
    %compute the gradient of the objective function with respect to the
    %choice variable z = [c_1(1) c_1(2) c_2(1)] using the gradients
    %computed above.  Note gradV is a 3x2 matrix as we have computed the
    %gradient for each of the two possible states
gradV=alpha(1).*psi.* c1.^(-sigma).*gradc1...
        +alpha(2).*psi.* c2.^(-sigma).*gradc2...
        -alpha(1).*(1-psi)./(1-l1).*gradl1...
        -alpha(2).*(1-psi)./(1-l2).*gradl2...
        +beta.*(V_x.*gradxprime+V_R.*gradRprime);
    
    
    %Sum over both states weighted by probabilities 
    grad =gradV*P(s_,:)';
    
    %A hack. Workers can't supply labor more than 1.  If this is the case
    %just return a large number and hope the root finder will move away
    %from this region
     if max([l1 l2]) >1
                grad=abs(z)+100;

     end
    %Sometimes we get complex numbers, again return a large number and hope
    %we move out of the region
    if ~isreal(grad)
    
    grad=abs(grad)+100;
    end
    else
        %if frac < 0 then return a large number and hope the rootfinder
        %moves out of this region.  Seems to work some of the time.
        grad=abs(z)+100;
   
    
    end
     

end