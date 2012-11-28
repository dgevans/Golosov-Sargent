function [l1 gradl1 l2 gradl2] = computeL(c1,gradc1,c2,gradc2,R,gradR,...
                                            theta_1,theta_2,g,n1,n2)
    
    %Get g in right format for computations
    g = kron(ones(3,1),g');
    %Compute l2 first
    l2 = ( n1*c1+n2*c2+g+n1*theta_2*R-n1*theta_1  )./(theta_2*(n2+R*n1));
    %now gradl2
    gradl2 = n1*gradR./(n2+n1*R) - n1*gradR.*l2./(n2+n1*R)...
             +n1*gradc1./(theta_2*(n2+n1*R)) +n2*gradc2./(theta_2*(n2+n1*R));
    %now l1
    l1 = 1 - (1-l2).*R*theta_2/theta_1;
    gradl1 = gradl2.*R*theta_2/theta_1 - (1-l2).*gradR*theta_2/theta_1;
         
         
end


