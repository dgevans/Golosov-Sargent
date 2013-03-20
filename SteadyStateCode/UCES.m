function [ u,uc,ul,ucc,ull ] = UCES( c , l , Para )
%UCES Summary of this function goes here
%   Detailed explanation goes here
    psi = Para.psi;
    u = psi*log(c)+(1-psi)*log(1-l);
    uc = psi./c;
    ucc = -psi*c.^(-2);
    ul = -(1-psi)./(1-l);
    ull = -(1-psi).*(1-l).^(-2);
end

