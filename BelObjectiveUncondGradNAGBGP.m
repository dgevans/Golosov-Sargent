function [ grad, user,iflag] = BelObjectiveUncondGradNAGBGP(n,x,user,iflag)
global V Vcoef R u2btild Par s_
%BELOBJECTIVEUNCOND Computes the Bellman objective with 
%   Detailed explanation goes here
    psi = Par.psi;
    beta =  Par.beta;
    P = Par.P;
    theta_1 = Par.theta(1);
    theta_2 = Par.theta(2);
    g = Par.g;
    alpha = Par.alpha;
    n1 = Par.n1;
    n2 = Par.n2;
    sigma = 1;
    
    frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
    
    if (min(x)>0 && frac>0)
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
    %compute components from unconstrained guess
    [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
%      h=.0001;
%  Der(1)=(computeC2_2(c1_1+h,c1_2,c2_1,R,s_,P,sigma)-computeC2_2(c1_1-h,c1_2,c2_1,R,s_,P,sigma))/(2*h)
%  Der(2)=(computeC2_2(c1_1,c1_2+h,c2_1,R,s_,P,sigma)-computeC2_2(c1_1,c1_2-h,c2_1,R,s_,P,sigma))/(2*h)
%  Der(3)=(computeC2_2(c1_1,c1_2,c2_1+h,R,s_,P,sigma)-computeC2_2(c1_1,c1_2,c2_1-h,R,s_,P,sigma))/(2*h)

    [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);
% computeLCheck = @(xx) computeL(xx(1),xx(2),xx(3),computeC2_2(xx(1),xx(2),xx(3),R,s_,P,sigma),grad_c2_2,theta_1,theta_2,g,n1,n2);
% computeLCheck([c1_1,c1_2 ,c2_1])
% Der(1,:)=(computeLCheck([c1_1+h,c1_2 ,c2_1])-computeLCheck([c1_1-h,c1_2 ,c2_1]))/(2*h)
% Der(2,:)=(computeLCheck([c1_1,c1_2+h ,c2_1])-computeLCheck([c1_1,c1_2-h ,c2_1]))/(2*h)
% Der(3,:)=(computeLCheck([c1_1,c1_2 ,c2_1+h])-computeLCheck([c1_1,c1_2 ,c2_1-h]))/(2*h)
    [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,psi,beta,P);
%  computeL1Check = @(xx) computeL1(xx(1),xx(2),xx(3),computeC2_2(xx(1),xx(2),xx(3),R,s_,P,sigma),grad_c2_2,theta_1,theta_2,g,n1,n2);
%   computeL1Check([c1_1,c1_2 ,c2_1])
%  computeL2Check = @(xx) computeL2(xx(1),xx(2),xx(3),computeC2_2(xx(1),xx(2),xx(3),R,s_,P,sigma),grad_c2_2,theta_1,theta_2,g,n1,n2);
%   computeL2Check([c1_1,c1_2 ,c2_1])
%  ComputeBtildeprimeCheck  =  @(xx) computeBtildeprime(xx(1),xx(2),xx(3),computeC2_2(xx(1),xx(2),xx(3),R,s_,P,sigma),grad_c2_2,computeL1Check([xx(1),xx(2),xx(3)]),computeL2Check([xx(1),xx(2),xx(3)]),l1grad,l2grad,u2btild,s_,psi,beta,P)
%  
%  ComputeBtildeprimeCheck([c1_1,c1_2 ,c2_1])
%  Der(1,:)=(ComputeBtildeprimeCheck([c1_1+h,c1_2 ,c2_1])-ComputeBtildeprimeCheck([c1_1-h,c1_2 ,c2_1]))/(2*h)
% Der(2,:)=(ComputeBtildeprimeCheck([c1_1,c1_2+h ,c2_1])-ComputeBtildeprimeCheck([c1_1,c1_2-h ,c2_1]))/(2*h)
%  Der(3,:)=(ComputeBtildeprimeCheck([c1_1,c1_2 ,c2_1+h])-ComputeBtildeprimeCheck([c1_1,c1_2 ,c2_1-h]))/(2*h)

 
 
 
    %compute objective
    grad1 = zeros(3,1);
    X = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
    Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi)...
            +beta*funeval(Vcoef{1},V(1),X));
    dV = funeval(Vcoef{1},V(1),X,eye(2));
      % Direct gardients with c_1_1,c_1_2,c_2_!
    grad1(1) = P(s_,1)*(alpha(1)*psi*c1_1^(-1)+beta*c2_1^(-1)*dV(2)); %<ok - Anmol>
    grad1(2) = 0; %<ok - Anmol>
    grad1(3) = P(s_,1)*(alpha(2)*psi*c2_1^(-1)-beta*c2_1^(-2)*(psi*btildprime(1)*dV(1)+c1_1*dV(2))); %<ok - Anmol>
    
    grad1 = grad1+P(s_,1)*( psi*grad_btildprime(:,1)*c2_1^(-1)*beta*dV(1) ...
        - alpha(1)*(1-psi)*l1grad(:,1)/(1-l1(1))...
        -alpha(2)*(1-psi)*l2grad(:,1)/(1-l2(1))); %<ok - Anmol>
    
    grad2 = zeros(3,1);
    X = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
    Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi)...
            +beta*funeval(Vcoef{2},V(2),X));
    dV = funeval(Vcoef{2},V(2),X,eye(2));
    grad2(1) = 0; %<ok - Anmol>
    grad2(2) = P(s_,2)*(alpha(1)*psi*c1_2^(-1)+beta*c2_2^(-1)*dV(2)); %<ok - Anmol>
    grad2(3) = 0;
    
    d_c2_2 = P(s_,2)*(alpha(2)*psi*c2_2^(-1)-beta*c2_2^(-2)*(dV(2)*c1_2 + psi*btildprime(2)*dV(1)));%<ok - Anmol>
    
    grad2 = grad2 + d_c2_2*grad_c2_2;
    grad2 = grad2 + P(s_,2)*( grad_btildprime(:,2)*psi*c2_2^(-1)*beta*dV(1)...
        -alpha(1)*(1-psi)*l1grad(:,2)/(1-l1(2))...
        -alpha(2)*(1-psi)*l2grad(:,2)/(1-l2(2)));
    
    grad = grad1+grad2;
     if max([l1 l2]) >1
                grad=abs(x)+100;

    end
    if ~isreal(grad)
    
    grad=abs(grad)+100;
    end
    else
        grad=abs(x)+100;
   
    
    end
     

end
% 
% 
% function [ c2_2 grad ] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma)
% 
%     %Compute c2_2 from formula
%     frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
%         /( P(s_,2) ); % <ok - Anmol>
%     c2_2 = frac^(-1/sigma); % <ok - Anmol>
%     grad=zeros(3,1);
%     %compute the gradients for c1_1,c1_2,c2_1
%     grad(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2)); % <ok - Anmol>
%     grad(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R; % <ok - Anmol>
%     grad(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2); % <ok - Anmol>
% end
% % CHECK GRADIENT
% 
% 
% 
% function [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
%     theta_1,theta_2,g,n1,n2)
% 
%     %Compute l1 form formula
%     l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1; % < ok - Anmol>
%     l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);  % < ok - Anmol>
%     l1(1) = l1_1num/l1_1den;  % < ok - Anmol>
%     l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2; % <ok - Anmol>
%     l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2); % <ok - Anmol>
%     l1(2) = l1_2num/l1_2den; % <ok - Anmol>
%     
%     %compute gradients of l1(1) for c1_1,c1_2,c2_1
%     l1grad(1,1) = (l1_1den*(n1-n2*theta_1*c2_1/c1_1^2)+l1_1num*n2*c2_1*theta_1/c1_1^2)/l1_1den^2; %<ok - Anmol>
%     l1grad(2,1) = 0;  % <ok - Anmol>
%     l1grad(3,1) = (l1_1den*(n2+n2*theta_1/c1_1)-l1_1num*n2*theta_1/c1_1)/l1_1den^2;  % <ok - Anmol>
%     
%     %compute gradients of l1(1) for c1_1,c1_2,c2_1
%     l1grad(1,2) = 0; % <ok - Anmol>
%     l1grad(2,2) = (l1_2den*(n1-n2*theta_1*c2_2/c1_2^2)+l1_2num*n2*c2_2*theta_1/c1_2^2)/l1_2den^2; % <ok - Anmol>
%     l1grad(3,2) = 0; % <ok - Anmol>
%     %use chain rule for c2_2
%     d_c2_2 = (l1_2den*(n2+n2*theta_1/c1_2)-l1_2num*n2*theta_1/c1_2)/l1_2den^2; % <ok - Anmol>
%     l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2; % <ok - Anmol>
%     
%     %compute l2 from formula
%     l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1; % <ok - Anmol>
%     l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
%     l2(1) = l2_1num/l2_1den; % <ok - Anmol>
%     l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2; % <ok - Anmol>
%     l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2; % <ok - Anmol>
%     l2(2) = l2_2num/l2_2den; % <ok - Anmol>
%     
%     %compute gradients of l2(1) for c1_1,c1_2,c2_1
%     l2grad(1,1) = (l2_1den*(n1+n1*theta_2/c2_1)-l2_1num*n1*theta_2/c2_1)/l2_1den^2;  % <ok - Anmol>
%     l2grad(2,1) = 0; % <ok - Anmol>
%     l2grad(3,1) = (l2_1den*(n2-n1*c1_1*theta_2/c2_1^2)+l2_1num*n1*c1_1*theta_2/c2_1^2)/l2_1den^2; % <ok - Anmol>
%     
%     %compute gradients of l2(2) for c1_1,c1_2,c2_1
%     l2grad(1,2) = 0; % <ok - Anmol>
%     l2grad(2,2) = (l2_2den*(n1+n1*theta_2/c2_2)-l2_2num*n1*theta_2/c2_2)/l2_2den^2; % <ok - Anmol>
%     l2grad(3,2) = 0; % <ok - Anmol>
%     %use chain rule to get the effect of c2_2
%     d_c2_2 = (l2_2den*(n2-n1*c1_2*theta_2/c2_2^2)+l2_2num*n1*c1_2*theta_2/c2_2^2)/l2_2den^2; % <ok - Anmol>
%     l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
%     
% end
% 
% function [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
%    u2btild,s_,psi,beta,P)
%     %get expected value of marginal utility of agent 2
%     Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
%  
%     %compute btildeprime from formula
%     btildprime(1) = u2btild/(beta*Eu2*psi)...
%         +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1))); % <Anmol - psi correction>
% 
%     %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
%     grad_btildprime(1,1) = 1-(1-psi)*l1(1)/(psi*(1-l1(1)));  % <ok - Anmol>
%     grad_btildprime(2,1) = 0;  % <ok - Anmol>
%     grad_btildprime(3,1) =u2btild*P(s_,1)*c2_1^(-2)/(beta*psi*Eu2^2)...  % <Anmol psi correction>
%         -1+(1-psi)*l2(1)/(psi*(1-l2(1)));  % <ok - Anmol>
% 
%     %figure out their affects through c2_2, l1_1,l2_1
%     d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(beta*psi*Eu2^2); % <Anmol psi correction>
%     d_l1_1 = -((1-psi)*c1_1/psi)/(1-l1(1))^2; % <ok - Anmol>
%     d_l2_1 = ((1-psi)*c2_1/psi)/(1-l2(1))^2;  % <ok - Anmol>
%     grad_btildprime(:,1) = grad_btildprime(:,1) + d_c2_2*grad_c2_2+d_l1_1*l1grad(:,1)+d_l2_1*l2grad(:,1); %<ok - Anmol>
% 
%     %Compute btildprime(2) from formula
%     btildprime(2) = u2btild/(psi*beta*Eu2)...
%         +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2))); %<Anmol psi correction>
% 
% 
%     %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
%     grad_btildprime(1,2) = 0; %<ok - Anmol>
%     grad_btildprime(2,2) = 1-(1-psi)*l1(2)/(psi*(1-l1(2))); %<ok - Anmol>
%     grad_btildprime(3,2) = u2btild*P(s_,1)*c2_1^(-2)/(psi*beta*Eu2^2);
%     %figure out their affects through c2_2, l1_2,l2_2
%     d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(psi*beta*Eu2^2)-1+(1-psi)*l2(2)/(psi*(1-l2(2)));
%     d_l1_2 = -((1-psi)*c1_2/psi)/(1-l1(2))^2;
%     d_l2_2 = ((1-psi)*c2_2/psi)/(1-l2(2))^2;
% 
%     grad_btildprime(:,2) = grad_btildprime(:,2) + d_c2_2*grad_c2_2+d_l1_2*l1grad(:,2)+d_l2_2*l2grad(:,2);
% 
% end
% 
% 
% function [l1] = computeL1(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
%     theta_1,theta_2,g,n1,n2)
% 
%     %Compute l1 form formula
%     l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1; % < ok - Anmol>
%     l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);  % < ok - Anmol>
%     l1(1) = l1_1num/l1_1den;  % < ok - Anmol>
%     l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2; % <ok - Anmol>
%     l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2); % <ok - Anmol>
%     l1(2) = l1_2num/l1_2den; % <ok - Anmol>
%     
%     %compute gradients of l1(1) for c1_1,c1_2,c2_1
%     l1grad(1,1) = (l1_1den*(n1-n2*theta_1*c2_1/c1_1^2)+l1_1num*n2*c2_1*theta_1/c1_1^2)/l1_1den^2; %<ok - Anmol>
%     l1grad(2,1) = 0;  % <ok - Anmol>
%     l1grad(3,1) = (l1_1den*(n2+n2*theta_1/c1_1)-l1_1num*n2*theta_1/c1_1)/l1_1den^2;  % <ok - Anmol>
%     
%     %compute gradients of l1(1) for c1_1,c1_2,c2_1
%     l1grad(1,2) = 0; % <ok - Anmol>
%     l1grad(2,2) = (l1_2den*(n1-n2*theta_1*c2_2/c1_2^2)+l1_2num*n2*c2_2*theta_1/c1_2^2)/l1_2den^2; % <ok - Anmol>
%     l1grad(3,2) = 0; % <ok - Anmol>
%     %use chain rule for c2_2
%     d_c2_2 = (l1_2den*(n2+n2*theta_1/c1_2)-l1_2num*n2*theta_1/c1_2)/l1_2den^2; % <ok - Anmol>
%     l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2; % <ok - Anmol>
%     
%     %compute l2 from formula
%     l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1; % <ok - Anmol>
%     l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
%     l2(1) = l2_1num/l2_1den; % <ok - Anmol>
%     l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2; % <ok - Anmol>
%     l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2; % <ok - Anmol>
%     l2(2) = l2_2num/l2_2den; % <ok - Anmol>
%     
%     %compute gradients of l2(1) for c1_1,c1_2,c2_1
%     l2grad(1,1) = (l2_1den*(n1+n1*theta_2/c2_1)-l2_1num*n1*theta_2/c2_1)/l2_1den^2;  % <ok - Anmol>
%     l2grad(2,1) = 0; % <ok - Anmol>
%     l2grad(3,1) = (l2_1den*(n2-n1*c1_1*theta_2/c2_1^2)+l2_1num*n1*c1_1*theta_2/c2_1^2)/l2_1den^2; % <ok - Anmol>
%     
%     %compute gradients of l2(2) for c1_1,c1_2,c2_1
%     l2grad(1,2) = 0; % <ok - Anmol>
%     l2grad(2,2) = (l2_2den*(n1+n1*theta_2/c2_2)-l2_2num*n1*theta_2/c2_2)/l2_2den^2; % <ok - Anmol>
%     l2grad(3,2) = 0; % <ok - Anmol>
%     %use chain rule to get the effect of c2_2
%     d_c2_2 = (l2_2den*(n2-n1*c1_2*theta_2/c2_2^2)+l2_2num*n1*c1_2*theta_2/c2_2^2)/l2_2den^2; % <ok - Anmol>
%     l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
%     
% end
% function [l2 ] = computeL2(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
%     theta_1,theta_2,g,n1,n2)
% 
%     %Compute l1 form formula
%     l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1; % < ok - Anmol>
%     l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);  % < ok - Anmol>
%     l1(1) = l1_1num/l1_1den;  % < ok - Anmol>
%     l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2; % <ok - Anmol>
%     l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2); % <ok - Anmol>
%     l1(2) = l1_2num/l1_2den; % <ok - Anmol>
%     
%     %compute gradients of l1(1) for c1_1,c1_2,c2_1
%     l1grad(1,1) = (l1_1den*(n1-n2*theta_1*c2_1/c1_1^2)+l1_1num*n2*c2_1*theta_1/c1_1^2)/l1_1den^2; %<ok - Anmol>
%     l1grad(2,1) = 0;  % <ok - Anmol>
%     l1grad(3,1) = (l1_1den*(n2+n2*theta_1/c1_1)-l1_1num*n2*theta_1/c1_1)/l1_1den^2;  % <ok - Anmol>
%     
%     %compute gradients of l1(1) for c1_1,c1_2,c2_1
%     l1grad(1,2) = 0; % <ok - Anmol>
%     l1grad(2,2) = (l1_2den*(n1-n2*theta_1*c2_2/c1_2^2)+l1_2num*n2*c2_2*theta_1/c1_2^2)/l1_2den^2; % <ok - Anmol>
%     l1grad(3,2) = 0; % <ok - Anmol>
%     %use chain rule for c2_2
%     d_c2_2 = (l1_2den*(n2+n2*theta_1/c1_2)-l1_2num*n2*theta_1/c1_2)/l1_2den^2; % <ok - Anmol>
%     l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2; % <ok - Anmol>
%     
%     %compute l2 from formula
%     l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1; % <ok - Anmol>
%     l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
%     l2(1) = l2_1num/l2_1den; % <ok - Anmol>
%     l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2; % <ok - Anmol>
%     l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2; % <ok - Anmol>
%     l2(2) = l2_2num/l2_2den; % <ok - Anmol>
%     
%     %compute gradients of l2(1) for c1_1,c1_2,c2_1
%     l2grad(1,1) = (l2_1den*(n1+n1*theta_2/c2_1)-l2_1num*n1*theta_2/c2_1)/l2_1den^2;  % <ok - Anmol>
%     l2grad(2,1) = 0; % <ok - Anmol>
%     l2grad(3,1) = (l2_1den*(n2-n1*c1_1*theta_2/c2_1^2)+l2_1num*n1*c1_1*theta_2/c2_1^2)/l2_1den^2; % <ok - Anmol>
%     
%     %compute gradients of l2(2) for c1_1,c1_2,c2_1
%     l2grad(1,2) = 0; % <ok - Anmol>
%     l2grad(2,2) = (l2_2den*(n1+n1*theta_2/c2_2)-l2_2num*n1*theta_2/c2_2)/l2_2den^2; % <ok - Anmol>
%     l2grad(3,2) = 0; % <ok - Anmol>
%     %use chain rule to get the effect of c2_2
%     d_c2_2 = (l2_2den*(n2-n1*c1_2*theta_2/c2_2^2)+l2_2num*n1*c1_2*theta_2/c2_2^2)/l2_2den^2; % <ok - Anmol>
%     l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
%     
% end