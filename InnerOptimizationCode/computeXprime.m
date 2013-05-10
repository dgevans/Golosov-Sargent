function [ xprime,gradxprime ] = computeXprime( c1,gradc1,c2,gradc2,Rprime,gradRprime,l1,gradl1,l2,gradl2,...
                                          P,sigma,psi,beta,s_,x)
%COMPUTEXPRIME %Computes the choice of the state variable xprime tomorrow in the
%standard 3x2 format as well as gradient with respect to z (note this
%is unfortunated notation, xprime is refering to xprime
%($u_{c,2}\tilde b'$), while x is the vector [c_1(1), c_1(2), c_2(1)]).
    

    %First create c2 alt.  Here c2alt is a matrix of the fomr
    %   c_2(2)  c_2(1)
    %   c_2(2)  c_2(1)
    %   c_2(2)  c_2(1)
    %This is so when we have multiplications like c2.*c2alt we get
    %   c_2(2)c_2(1) c_2(1)c_2(2)
    %   c_2(2)c_2(1) c_2(1)c_2(2)
    %   c_2(2)c_2(1) c_2(1)c_2(2)
    S = length(P(1,:));
    %Now the expected marginal utility of agent 2.  Again want it in 3x2
    %format
    Euc2 = kron(ones(1,S),psi*c2.^(-sigma)*(P(s_,:)'));
    
 
    Palt = kron(ones(2*S-1,1),P(s_,:));
    %gives a 2S-1 x S matrix with each column representing the partial derivative of
    %Euc2 with eash S
    grad_Euc2 = sigma*psi*c2.^(-sigma-1).*Palt.*gradc2;
    grad_Euc2 = sum(grad_Euc2,2);
    %copy this vector for each S
    grad_Euc2 = kron(grad_Euc2,ones(1,S));
    
    %create new 3x2 P and Palt
    %P = kron(ones(2*S-1,1),P(s_,:));
    %Palt = fliplr(P);
    % x*psi*c2.^(-sigma)./(beta*Euc2) + (1-psi)*l2./(1-l2)...
    %         -(1-psi)*Rprime.*l1./(1-l1)+psi*c1.*c2.^(-sigma)-psi*c2.^(1-sigma);
   
    %Now compute xprime from formula in notes
    xprime_beta = x*psi*c2.^(-sigma)./(Euc2) + (1-psi)*l2./(1-l2)...
             -(1-psi)*Rprime.*l1./(1-l1)+psi*c1.*c2.^(-sigma)-psi*c2.^(1-sigma);
    xprime=xprime_beta./beta;
         %Now compute the gradient
    gradxprime = ( -sigma*x*psi*c2.^(-sigma-1)./(Euc2)...              % x*psi*c2.^(-sigma)./(beta*Euc2) with c2
                   -sigma*psi*c2.^(-sigma-1).*c1...                         % psi*c1.*c2.^(-sigma) with c2
                   -(1-sigma)*psi*c2.^(-sigma)).*gradc2...                  % psi*c2.^(1-sigma) with c2
                +(x*psi*c2.^(-sigma)./((Euc2).^2) ).*grad_Euc2... % x*psi*c2.^(-sigma)./(beta*Euc2) with Euc2
                +psi*c2.^(-sigma).*gradc1...                                %  psi*c1.*c2.^(-sigma) with c1
                +(1-psi)*gradl2./((1-l2).^2)...                             % (1-psi)*l2./(1-l2) with l2
                -(1-psi)*Rprime.*gradl1./((1-l1).^2)...                     % -(1-psi)*Rprime.*l1./(1-l1) with l1
                -(1-psi)*l1.*gradRprime./(1-l1);                            % -(1-psi)*Rprime.*l1./(1-l1) with Rprime
gradxprime=gradxprime./beta;
end

