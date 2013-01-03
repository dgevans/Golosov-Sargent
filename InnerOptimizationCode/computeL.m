function [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2)
%COMPUTEL computes l_1 and l_2, the labor supply  of agent 1 and 2 in the
%standard 3x2 format, along with their gradients with respect to z.  Uses
%c1, c2, Rprime computed using computeC2 and computeRprime as well as their
%gradients.  Also passed are the primitives theta_1,theta_2, n_1, n_2 and
%the vector of government expenditures g.
    
    %Get g in right format for computations
    g = g(:)';
    S = max([length(g),length(theta_1),length(theta_2)]);
    if(length(g) > 1)
        g = kron(ones(2*S-1,1),g);
    end
    if(length(theta_1) > 1)
        theta_1 = theta_1(:)';
        theta_1 = kron(ones(2*S-1,1),theta_1);
    end
    if(length(theta_2) > 1)
        theta_2 = theta_2(:)';
        theta_2 = kron(ones(2*S-1,1),theta_2);
    end
    %Compute l2 first
    l2 = ( n1*c1+n2*c2+g+n1.*theta_2*Rprime-n1.*theta_1  )./(theta_2.*(n2+Rprime*n1));
    %now gradl2
    gradl2 = n1*gradRprime./(n2+n1*Rprime) - n1*gradRprime.*l2./(n2+n1*Rprime)...
             +n1*gradc1./(theta_2.*(n2+n1*Rprime)) +n2*gradc2./(theta_2.*(n2+n1*Rprime));
    %now l1
    l1 = 1 - (1-l2).*Rprime.*theta_2./theta_1;
    gradl1 = gradl2.*Rprime.*theta_2./theta_1 - (1-l2).*gradRprime.*theta_2./theta_1;
         
         
end


