function [ res, iflag] = resFOCBGP_alt(n,x,iflag)
global V Vcoef R u2btild Par s_ flagCons
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
    u2btildLL=Par.u2btildLL;
    u2btildUL=Par.u2btildUL;
    sigma = 1;
    frac = (R*P(s_,1)*x(1)^(-1)+R*P(s_,2)*x(2)^(-1)-P(s_,1)*x(3)^(-1))...
        /( P(s_,2) );
    
    
    
    if (min(x(1:3))>0 && frac>0)
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
     MuL = zeros(2,1);
     MuH = zeros(2,1);
    switch flagCons
    case 'LL_'
       % lower limit binds for state 1 only
       MuL(1)=x(4);
       MuL(2)=0;
       u2btildprime(1)=u2btildLL;
       u2btildprime(2)=x(5);
       
    case '_LL'
       % lower limit binds for state 2 only
       MuL(1)=0;
       MuL(2)=x(5);
       u2btildprime(1)=x(4);
       u2btildprime(2)=u2btildLL;
       
    case 'LLLL'
      % lower limit binds for both the states
       MuL(1)=x(4);
       MuL(2)=x(5);
       u2btildprime(1)=u2btildLL;   
       u2btildprime(2)=u2btildLL;     
        
    case 'UL_'
     % upper limit binds for state 1 only

       MuH(1)=x(4);
       MuH(2)=0;
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=x(5);
       
        
    case '_UL'
         % upper limit binds for state 2 only
       MuH(1)=0;
       MuH(2)=x(5);
       u2btildprime(1)=x(4);
       u2btildprime(2)=u2btildUL;
        
        
    case 'ULUL'
        
        
       % upper limit binds for both the states
       MuH(1)=x(4);
       MuH(2)=x(5);
       u2btildprime(1)=u2btildUL;
       u2btildprime(2)=u2btildUL;     
        
    otherwise
       MuL(1)=0;
       MuL(2)=0;
       u2btildprime(1)=x(4);
       u2btildprime(2)=x(5);     
        
    end
    lambda_I(1) = x(6);
    lambda_I(2) = x(7);
    res = zeros(7,1);
    
    %compute components from unconstrained guess
    [c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
    
    [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);
    
[btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,psi,beta,P);
    
    %compute objective
    grad1 = zeros(3,1);
    X = [u2btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
    Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi)...
            +beta*funeval(Vcoef{1},V(1),X));
    % derivatives of the value function with x and R for s=1
    dV = funeval(Vcoef{1},V(1),X,eye(2));
    
    
    % Direct gardients with c_1_1,c_1_2,c_2_!
    grad1(1) = P(s_,1)*(alpha(1)*psi*c1_1^(-1)+beta*c2_1^(-1)*dV(2));
    grad1(2) = 0;
    grad1(3) = P(s_,1)*(alpha(2)*psi*c2_1^(-1) ...
        -c2_1^(-2)*(psi*btildprime(1)*(-lambda_I(1)) + c1_1*beta*dV(2))); % Anmol shifted the beta to dV(2)
% Indirect gradients via btildprime, l1 and l2

    grad1 = grad1+P(s_,1)*( psi*grad_btildprime(:,1)*c2_1^(-1)*(-lambda_I(1))... % ok
        - alpha(1)*(1-psi)*l1grad(:,1)/(1-l1(1))... %ok
        -alpha(2)*(1-psi)*l2grad(:,1)/(1-l2(1))); %ok
   % FOC with respect to x'(1)
    res(4) = P(s_,1)*lambda_I(1)+P(s_,1)*beta*dV(1)+MuL(1)-MuH(1); % Anmol - Changed the sign on MuH to minus
   %ok lambda(1)=dV(1)
   
    % Definition of x'
    res(6) = u2btildprime(1)-psi*c2_1^(-1)*btildprime(1);

    grad2 = zeros(3,1);
    X = [u2btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
    Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi)...
            +beta*funeval(Vcoef{2},V(2),X));
    dV = funeval(Vcoef{2},V(2),X,eye(2));
    % Direct gradients for s=2
    grad2(1) = 0; % c_1_1
    grad2(2) = P(s_,2)*(alpha(1)*psi*c1_2^(-1)+beta*c2_2^(-1)*dV(2));  % c_1_2
    grad2(3) = 0; % c_2_1
    
    % Now take derivatives with c_2_2
    d_c2_2 = P(s_,2)*(alpha(2)*psi*c2_2^(-1)-c2_2^(-2)*(beta*dV(2)*c1_2+psi*btildprime(2)*(-lambda_I(2))));
    grad2 = grad2 + d_c2_2*grad_c2_2;

    % Indirect gradients via btildprime, l1 and l2    
    grad2 = grad2 + P(s_,2)*( grad_btildprime(:,2)*psi*c2_2^(-1)*(-lambda_I(2))...
        -alpha(1)*(1-psi)*l1grad(:,2)/(1-l1(2))...
        -alpha(2)*(1-psi)*l2grad(:,2)/(1-l2(2)));
    % FOC with x'(2)
    res(5) = P(s_,2)*lambda_I(2)+P(s_,2)*beta*dV(1)+MuL(2)-MuH(2); % Anmol - Changed the sign on MuH to minus and corrected typo for dV(2-->1)
    % Definition of x'(2)
    res(7) = u2btildprime(2)-psi*c2_2^(-1)*btildprime(2);
    res(1:3) = grad1+grad2;
           
     if max([l1 l2]) >1
                res=abs(x)+100;

    end
     
    if ~isreal(grad1+grad2)
    
    res=abs(x)+100;
    end
    else
        res=abs(x)+100;
    end

end
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
% % 
% % 
% % function [ c2_2 grad ] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma)
% % 
% %     %Compute c2_2 from formula
% %     frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
% %         /( P(s_,2) );
% %     c2_2 = frac^(-1/sigma);
% %     grad=zeros(3,1);
% %     %compute the gradients for c1_1,c1_2,c2_1
% %     grad(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2));
% %     grad(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R;
% %     grad(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2);
% % end
% % 
% % function [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
% %     theta_1,theta_2,g,n1,n2)
% % 
% %     %Compute l1 form formula
% %     l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1;
% %     l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);
% %     l1(1) = l1_1num/l1_1den;
% %     l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2;
% %     l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2);
% %     l1(2) = l1_2num/l1_2den;
% %     
% %     %compute gradients of l1(1) for c1_1,c1_2,c2_1
% %     l1grad(1,1) = (l1_1den*(n1-n2*theta_1*c2_1/c1_1^2)+l1_1num*n2*c2_1*theta_1/c1_1^2)/l1_1den^2;
% %     l1grad(2,1) = 0;
% %     l1grad(3,1) = (l1_1den*(n2+n2*theta_1/c1_1)-l1_1num*n2*theta_1/c1_1)/l1_1den^2;
% %     
% %     %compute gradients of l1(2) for c1_1,c1_2,c2_1
% %     l1grad(1,2) = 0;
% %     l1grad(2,2) = (l1_2den*(n1-n2*theta_1*c2_2/c1_2^2)+l1_2num*n2*c2_2*theta_1/c1_2^2)/l1_2den^2;
% %     l1grad(3,2) = 0;
% %     %use chain rule for c2_2
% %     d_c2_2 = (l1_2den*(n2+n2*theta_1/c1_2)-l1_2num*n2*theta_1/c1_2)/l1_2den^2;
% %     l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2;
% %     
% %     %compute l2 from formula
% %     l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1;
% %     l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
% %     l2(1) = l2_1num/l2_1den;
% %     l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2;
% %     l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2;
% %     l2(2) = l2_2num/l2_2den;
% %     
% %     %compute gradients of l2(1) for c1_1,c1_2,c2_1
% %     l2grad(1,1) = (l2_1den*(n1+n1*theta_2/c2_1)-l2_1num*n1*theta_2/c2_1)/l2_1den^2;
% %     l2grad(2,1) = 0;
% %     l2grad(3,1) = (l2_1den*(n2-n1*c1_1*theta_2/c2_1^2)+l2_1num*n1*c1_1*theta_2/c2_1^2)/l2_1den^2;
% %     
% %     %compute gradients of l2(2) for c1_1,c1_2,c2_1
% %     l2grad(1,2) = 0;
% %     l2grad(2,2) = (l2_2den*(n1+n1*theta_2/c2_2)-l2_2num*n1*theta_2/c2_2)/l2_2den^2;
% %     l2grad(3,2) = 0;
% %     %use chain rule to get the effect of c2_2
% %     d_c2_2 = (l2_2den*(n2-n1*c1_2*theta_2/c2_2^2)+l2_2num*n1*c1_2*theta_2/c2_2^2)/l2_2den^2;
% %     l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
% %     
% % end
% % 
% % function [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
% %    u2btild,s_,psi,beta,P)
% %     %get expected value of marginal utility of agent 2
% %     Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
% %     
% %     %compute btildeprime from formula
% %     btildprime(1) = u2btild/(psi*beta*Eu2)...
% %         +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1)));
% % 
% %     %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
% %     grad_btildprime(1,1) = 1-(1-psi)*l1(1)/(psi*(1-l1(1)));
% %     grad_btildprime(2,1) = 0;
% %     grad_btildprime(3,1) =u2btild*P(s_,1)*c2_1^(-2)/(psi*beta*Eu2^2)...
% %         -1+(1-psi)*l2(1)/(psi*(1-l2(1)));
% % 
% %     %figure out their affects through c2_2, l1_1,l2_1
% %     d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(psi*beta*Eu2^2);
% %     d_l1_1 = -((1-psi)*c1_1/psi)/(1-l1(1))^2;
% %     d_l2_1 = ((1-psi)*c2_1/psi)/(1-l2(1))^2;
% %     grad_btildprime(:,1) = grad_btildprime(:,1) + d_c2_2*grad_c2_2+d_l1_1*l1grad(:,1)+d_l2_1*l2grad(:,1);
% % 
% %     %Compute btildprime(2) from formula
% %     btildprime(2) = u2btild/(psi*beta*Eu2)...
% %         +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2)));
% % 
% % 
% %     %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
% %     grad_btildprime(1,2) = 0;
% %     grad_btildprime(2,2) = 1-(1-psi)*l1(2)/(psi*(1-l1(2)));
% %     grad_btildprime(3,2) = u2btild*P(s_,1)*c2_1^(-2)/(beta*Eu2^2*psi);
% %     %figure out their affects through c2_2, l1_2,l2_2
% %     d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(beta*Eu2^2*psi)-1+(1-psi)*l2(2)/(psi*(1-l2(2)));
% %     d_l1_2 = -((1-psi)*c1_2/psi)/(1-l1(2))^2;
% %     d_l2_2 = ((1-psi)*c2_2/psi)/(1-l2(2))^2;
% % 
% %     grad_btildprime(:,2) = grad_btildprime(:,2) + d_c2_2*grad_c2_2+d_l1_2*l1grad(:,2)+d_l2_2*l2grad(:,2);
% % 
% % end
