function [ u,uc,ul,ucc,ull ] = UMix( c, l, Para )
%UMIX Summary of this function goes here
%   Detailed explanation goes here
    sigma = Para.sigma;
    psi = Para.psi;
    if sigma ==1
        u = psi*log(c) + (1-psi)*log(1-l);
    else
        u = psi*(c.^(1-sigma)-1)/(1-sigma) + (1-psi)*log(1-l); 
    end
    
    uc = psi* c.^(-sigma);
    ucc = -sigma*psi*c.^(-sigma-1);
    ul = -(1-psi)./(1-l);
    ull = - (1-psi).*(1-l).^(-2);

end

