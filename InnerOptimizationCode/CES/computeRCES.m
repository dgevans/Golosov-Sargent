function [ Rprime,gradRprime ] = computeRCES( c1,c2,gradc1,gradc2,sigma)
%COMPUTER Computes Rprime c_2(s)^(-sigma)/c_1(s)^(-sigma) in the 3x2 matrix
%format as well as its gradient with respect to z.  Uses the variables c1, c2
%computed in computeC2 as well as their gradients.  Also needed is the
%primitive sigma.

    Rprime = (c2.^(-sigma) )./(c1.^(-sigma));
    
    gradRprime = sigma*c2.^(-sigma).*c1.^(sigma-1).*gradc1...
           -sigma*c2.^(-sigma-1).*c1.^(sigma).*gradc2;
end

