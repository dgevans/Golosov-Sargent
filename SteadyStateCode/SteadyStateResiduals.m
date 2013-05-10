% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [res c1 c2 l1 l2]=SteadyStateResiduals(z,u2bdiff,RR,Para,s_)
Par=Para;
u2btild=u2bdiff;
R=RR;
n1=Para.n1;
n2=Para.n2;

%% GET THE Policy Rules
psi= Par.psi;beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta_1;
theta_2 = Par.theta_2;
g = Par.g;
sigma = Par.sigma;

    z = z(:)';

    S = length(P(1,:));
    c1 = z(1:S);
    c2_ = z(S+1:2*S-1);
    P_ = P(s_,:); P_(S) = [];
    frac = (R*dot(P(s_,:),c1.^(-sigma)) - dot(P_,c2_.^(-sigma)))/P(s_,S);
    %it is necessary for frac to be positive in order for the constraints
    %to be satisfied.
  if numel(beta)==1
    beta=ones(2*S-1,S)*beta;
    else
    beta=repmat(beta,2*S-1,1);
    end
if (min(z)>0 && frac>0)

    %compute components from unconstrained guess
    [c1,c2,gradc1,gradc2] = computeC2_2(c1,c2_,R,s_,P,sigma);
    [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                                theta_1,theta_2,g,n1,n2);
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                              P,sigma,psi,beta,s_,u2btild);

    % State next period
    xprime = xprime(1,:);
    Rprime = Rprime(1,:);


    res(1:S) = xprime-u2bdiff;
    res(S+1:2*S-1) = Rprime(1:S-1)-R;

    c1 = c1(1,:);
    c2 = c2(1,:);
    l1 = l1(1,:);
    l2 = l2(1,:);
        if max([l1 l2]) >1
                    res=abs(z)+100;

        end
        if ~isreal(res)

        res=abs(res)+100;
        end
        else
            res=abs(z)+100;


end


end
%
%res(4)=X(1,2)-R;

