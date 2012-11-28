function [ u ] = uAlt( c,l,psi,sigma )
%UALT Summary of this function goes here
%   Detailed explanation goes here
    if(sigma == 1)
        u = psi*log(c)+(1-psi)*log(1-l);
    else
        u = psi*c.^(1-sigma)/(1-sigma)+(1-psi)*log(1-l);
        

end

