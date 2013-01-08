% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [res c1 c2 l1 l2]=SteadyStateResiduals(x,u2bdiff,RR,Para,s)
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
s_=s;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

%% GET THE Policy Rules
psi= Par.psi;beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;
sigma = Par.sigma;

   frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );

if (min(x)>0 && frac>0)

    c1_1=x(1);
    c1_2=x(2);
    c2_1=x(3);

    %compute components from unconstrained guess
    [c1,c2,gradc1,gradc2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma);
    [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                                theta_1,theta_2,g,n1,n2);
    [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                              P,sigma,psi,beta,s_,u2btild);

    % State next period
    xprime = xprime(1,:);
    Rprime = Rprime(1,:);


    res(1)=xprime(1) -xprime(2);
    res(2)=Rprime(1) - Rprime(2);
    res(3)=xprime(1)-u2btild;

    c1 = c1(1,:);
    c2 = c2(1,:);
    l1 = l1(1,:);
    l2 = l2(1,:);
    if max([l1 l2]) >1
        res=abs(x)+100;
    end

    if ~isreal(res)
        res=abs(res)+100;
    end

else
    res=abs(x)+100;


end


end
%
%res(4)=X(1,2)-R;

