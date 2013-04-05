function [ u ] = uAltCES( c,l,gamma,sigma )
%UALT Utility function allowing for different values of risk adversion
    if(sigma == 1)
        u = log(c)-l.^(1+gamma)/(1+gamma);
    else
        u = c.^(1-sigma)/(1-sigma)-l.^(1+gamma)/(1+gamma);
    end
        

end

