function [ u,uc,ul,ucc,ull ] = UCRRA( c,l,Para )
%UCRRA Summary of this function goes here
%   Detailed explanation goes here
    sigma = Para.sigma;
    gamma = Para.gamma;
    u = c.^(1-sigma)./(1-sigma) - l.^(gamma+1)./(gamma+1);
    
    uc = c.^(-sigma);
    ucc = -sigma* c.^(-sigma-1);
    
    ul = -l.^gamma;
    ull = -gamma*l.^(gamma-1);

end

