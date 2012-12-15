function [ Rprime,gradRprime ] = computeR( c1,c2,gradc1,gradc2,sigma)
%COMPUTER Computes R and it's gradient for unconstrained optimization
%   Detailed explanation goes here

    Rprime = (c2.^(-sigma) )./(c1.^(-sigma));
    
    gradRprime = sigma*c2.^(-sigma).*c1.^(sigma-1).*gradc1...
           -sigma*c2.^(-sigma-1).*c1.^(sigma).*gradc2;
end

