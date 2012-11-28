function [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,Rprime,gradRprime,...
                                            theta_1,theta_2,g,n1,n2)
    
    %Get g in right format for computations
    if(size(g,1) > 1)
        g = g';
    end
    g = kron(ones(3,1),g);
    %Compute l2 first
    l2 = ( n1*c1+n2*c2+g+n1*theta_2*Rprime-n1*theta_1  )./(theta_2*(n2+Rprime*n1));
    %now gradl2
    gradl2 = n1*gradRprime./(n2+n1*Rprime) - n1*gradRprime.*l2./(n2+n1*Rprime)...
             +n1*gradc1./(theta_2*(n2+n1*Rprime)) +n2*gradc2./(theta_2*(n2+n1*Rprime));
    %now l1
    l1 = 1 - (1-l2).*Rprime*theta_2/theta_1;
    gradl1 = gradl2.*Rprime*theta_2/theta_1 - (1-l2).*gradRprime*theta_2/theta_1;
         
         
end


