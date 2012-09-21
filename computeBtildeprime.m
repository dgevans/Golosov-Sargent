function [btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
   u2btild,s_,psi,beta,P)
    %get expected value of marginal utility of agent 2
    Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
 
    %compute btildeprime from formula
    btildprime(1) = u2btild/(beta*Eu2*psi)...
        +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1))); % <Anmol - psi correction>

    %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
    grad_btildprime(1,1) = 1-(1-psi)*l1(1)/(psi*(1-l1(1)));  % <ok - Anmol>
    grad_btildprime(2,1) = 0;  % <ok - Anmol>
    grad_btildprime(3,1) =u2btild*P(s_,1)*c2_1^(-2)/(beta*psi*Eu2^2)...  % <Anmol psi correction>
        -1+(1-psi)*l2(1)/(psi*(1-l2(1)));  % <ok - Anmol>

    %figure out their affects through c2_2, l1_1,l2_1
    d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(beta*psi*Eu2^2); % <Anmol psi correction>
    d_l1_1 = -((1-psi)*c1_1/psi)/(1-l1(1))^2; % <ok - Anmol>
    d_l2_1 = ((1-psi)*c2_1/psi)/(1-l2(1))^2;  % <ok - Anmol>
    grad_btildprime(:,1) = grad_btildprime(:,1) + d_c2_2*grad_c2_2+d_l1_1*l1grad(:,1)+d_l2_1*l2grad(:,1); %<ok - Anmol>

    %Compute btildprime(2) from formula
    btildprime(2) = u2btild/(psi*beta*Eu2)...
        +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2))); %<Anmol psi correction>


    %compute grad of btildprime(1) with respect to c1_1,c1_2,c2_1
    grad_btildprime(1,2) = 0; %<ok - Anmol>
    grad_btildprime(2,2) = 1-(1-psi)*l1(2)/(psi*(1-l1(2))); %<ok - Anmol>
    grad_btildprime(3,2) = u2btild*P(s_,1)*c2_1^(-2)/(psi*beta*Eu2^2);
    %figure out their affects through c2_2, l1_2,l2_2
    d_c2_2 = u2btild*P(s_,2)*c2_2^(-2)/(psi*beta*Eu2^2)-1+(1-psi)*l2(2)/(psi*(1-l2(2)));
    d_l1_2 = -((1-psi)*c1_2/psi)/(1-l1(2))^2;
    d_l2_2 = ((1-psi)*c2_2/psi)/(1-l2(2))^2;

    grad_btildprime(:,2) = grad_btildprime(:,2) + d_c2_2*grad_c2_2+d_l1_2*l1grad(:,2)+d_l2_2*l2grad(:,2);

end