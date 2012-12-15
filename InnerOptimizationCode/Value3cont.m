function [ minusVobj,minusGrad] = Value3cont(z)
global V Vcoef R u2btild Par s_
%BELOBJECTIVEUNCOND Computes minus one times the gradient and value of the
%bellman equation objective with
%respect to c_1(1), c_1(2) and c_2(1).  Substitues out for the rest of
%the variables using their respective gradients.  z is the vector
%containing these three variables
%  Note that we pass global variables as well (because NAG did not allow us
%  to pass user defined information) V is a struct containing the
%  interpolation information.  Vcoef has coefficents for the value
%  function. R, u2btild and s_ are the previous period states.  Par is a
%  struct containig all the relevent parameters.  We return the gradient
%  grad.  user and iflag our variables that nag requires but we don't use.
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
    
    frac = (R*P(s_,1)*z(1)^(-sigma)+R*P(s_,2)*z(2)^(-sigma)-P(s_,1)*z(3)^(-sigma))...
        /( P(s_,2) );
        %it is necessary for frac to be positive in order for the constraints
    %to be satisfied.
    if (min(z)>0 && frac>0)
     c1_1=z(1);
     c1_2=z(2);
     c2_1=z(3);
    %compute components from unconstrained guess
    
    
    [c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
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
                                          P,sigma,psi,beta,s_,u2btild);
    %Computes the choice of the state variable xprime tomorrow in the
    %standard 3x2 format as well as gradient with respect to z (note this
    %is unfortunated notation, xprime is refering to u2btildprime
    %($u_{c,2}\tilde b'$), while z is the vector [c_1(1), c_1(2), c_2(1)].
    %Sorry for the confusion
 
    %compute objective and also the gradient of the object.  So first
    %compute continuation values and their gradient with respect to the
    %state variables
    Vprime(:,1) = funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)])*ones(3,1);
    Vprime(:,2) = funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)])*ones(3,1);
    V_x(:,1)=funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)],[1,0])*ones(3,1);
    V_x(:,2)=funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)],[1,0])*ones(3,1);
    V_R(:,1)=funeval(Vcoef{1},V(1),[xprime(1,1) Rprime(1,1)],[0,1])*ones(3,1);
    V_R(:,2)=funeval(Vcoef{2},V(2),[xprime(1,2) Rprime(1,2)],[0,1])*ones(3,1);
    
    %This is the value of the objective function for both states.  Note
    %this is a 3x2 matrix.  We will only use the first row later as all the
    %rows are the same
    Vrhs = alpha(1)*uAlt(c1,l1,psi,sigma)+alpha(2)*uAlt(c2,l2,psi,sigma) + beta*Vprime;
    

    
    %compute the gradient of the objective function with respect to the
    %choice variable z = [c_1(1) c_1(2) c_2(1)] using the gradients
    %computed above.  Note gradV is a 3x2 matrix as we have computed the
    %gradient for each of the two possible states
gradV=alpha(1).*psi.* c1.^(-sigma).*gradc1...
        +alpha(2).*psi.* c2.^(-sigma).*gradc2...
        -alpha(1).*(1-psi)./(1-l1).*gradl1...
        -alpha(2).*(1-psi)./(1-l2).*gradl2...
        +beta*(V_x.*gradxprime+V_R.*gradRprime);
    
    %Compute expected value using the transition matrix
    minusGrad =-gradV*P(s_,:)';
    minusVobj = -Vrhs(1,:)*P(s_,:)';
    
    %code to return large values if in a bad region
     if max([l1 l2]) >1
                grad=abs(z)+100;

    end
    if ~isreal(minusGrad)
    minusVobj = 100;
    minusGrad=-abs(minusGrad)-100;
    end
    else
        minusVobj = 100;
        minusGrad=-abs(z)-100;
   
    
    end
     

end