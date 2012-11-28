function [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,u2btild)
%COMPUTEXPRIME Computes xprime and its gradient.
%   It is assumed that all the passed variables except for P are either
%   scalars or have dimension 3x2
    

    %First create c2 alt 
    c2alt = fliplr(c2);
    gradc2alt = fliplr(gradc2);
    %Now Ec2 remember want it to be 3x2
    Euc2 = kron(ones(1,2),c2.^(-sigma)*(P(s_,:)'));
    
    %create new 3x2 P and Palt
    P = kron(ones(3,1),P(s_,:));
    Palt = fliplr(P);
    
    %Now compute xprime
    xprime = u2btild*psi*c2.^(-sigma)./(beta*Euc2) + (1-psi)*l2./(1-l2)...
             -(1-psi)*Rprime.*l1./(1-l1)+psi*c1.*c2.^(-sigma)-psi*c2.^(1-sigma);
    %Now compute the gradient
    gradxprime = ( -sigma*u2btild*psi*c2.^(-sigma-1)./(beta*Euc2)...
                    + (sigma*u2btild*psi*c2.^(-2*sigma-1).*P.*beta)./((beta*Euc2).^2)...
                   -sigma*psi*c2.^(-sigma-1).*c1-(1-sigma)*psi*c2.^(-sigma)).*gradc2...
                +(sigma*u2btild*psi*c2.^(-sigma).*c2alt.^(-sigma-1).*beta.*Palt)...
                 ./((beta*Euc2).^2)*gradc2alt+psi*c2^(-sigma).*gradc1...
                +(1-psi)*gradl2./((1-l2).^2)-(1-psi)*Rprime.*gradl1./((1-l1).^2)...
                -(1-psi)*l1.*gradRprime./(1-l1);
end

