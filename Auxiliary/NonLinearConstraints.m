function [cineq,ceq,grad_ineq grad_eq] = NonLinearConstraints(x,u2btild,R,s_,Par)
grad_ineq=[];
cineq=[];
theta_1=Par.theta_1;
sigma=Par.sigma;
gamma=Par.gamma;
g=Par.g;
beta=Par.beta;
P=Par.P;
c1(1)=x(1); %consumption of agent 1 state 1
c1(2)=x(2); %consumption of agent 1 state 2
c2(1)=x(3);  %consumption of agent 2 state 1
c2(2)=x(4); %consumption of agent 2 state 2
l1(1)=x(5); %labor supply of agent 1 state 1
l1(2)=x(6); %labor supply of agent 1 state 2
u2btildprime(1)=x(7);
u2btildprime(2)=x(8);

    % Get the expected value of the marginal utilities 
    Eu1=P(s_,1)*c1(1)^(-sigma)+P(s_,2)*c1(2)^(-sigma);
        Eu2=P(s_,1)*c2(1)^(-sigma)+P(s_,2)*c2(2)^(-sigma);

ceq=[(c2(1)-c1(1))+l1(1)^(1+gamma)/c1(1)^(-sigma)+u2btildprime(1)/c2(1)^(-sigma)-u2btild/(beta*Eu2);...         % Implementability state 1
    (c2(2)-c1(2))+l1(2)^(1+gamma)/c1(2)^(-sigma)+u2btildprime(2)/c2(2)^(-sigma)-u2btild/(beta*Eu2);...             % Implementability state 2
Eu2/Eu1-R; ...                                                                                   % Bond Pricing   
c1(1)+c2(1)+g(1)-theta_1*l1(1); ...                                                               % Resource Constraint state 1
c1(2)+c2(2)+g(2)-theta_1*l1(2);]; 


% Derivative of the Implementability constraint
    dI(:,1)=[-1+l1(1)^(1+gamma)*sigma*c1(1)^(sigma-1); ...                                                               % c1(1)
         0; ...                                                                                                          % c1(2)
        1+u2btildprime(1)*sigma*c2(1)^(sigma-1)-((sigma*u2btild*P(s_,1)*c2(1)^(-sigma-1))/(beta*Eu2^2)); ....             %c2(1)
        -(u2btild*P(s_,2)*sigma*c2(2)^(-sigma-1)/(beta*Eu2^2));...                                                         %c2(2)
         (1+gamma)*l1(1)^gamma*c1(1)^(sigma);...                                                                         %l1(1)
         0; ...                                                                                                           %l1(2)
         c2(1)^(sigma); ...                                                                                             %u2btildprime(1)
         0;];                                                                                                           %u2btldprime (2)   
        
     
    dI(:,2)=[0; ...                                                                                                     %c1(1)
        -1+l1(2)^(1+gamma)*sigma*c1(2)^(sigma-1); ...                                                                   %c1(2)
        -(u2btild*P(s_,1)*sigma*c2(1)^(-sigma-1)/(beta*Eu2^2));   ...                                                     %c2(1)
        1+u2btildprime(2)*sigma*c2(2)^(sigma-1)-((sigma*u2btild*P(s_,2)*c2(2)^(-sigma-1))/(beta*Eu2^2)); ....            %c2(2)
        0;                                                                                                              %l1(1)
       (1+gamma)*l1(2)^gamma*c1(2)^(sigma);...                                                                          %l1(2)
          0;                                                                                                            %u2btildprime(1)
         c2(2)^(sigma);];                                                                                            %u2btildprime(2)
         
     
     % Derivative of the Bond Pricing constraint
 dB=[ (Eu2/Eu1^2)*P(s_,1)*sigma*c1(1)^(-sigma-1);...                                            %c1(1)
      (Eu2/Eu1^2)*P(s_,2)*sigma*c1(2)^(-sigma-1);...                                            %c1(2)
      (-1/Eu1)*P(s_,1)*sigma*c2(1)^(-sigma-1);...                                               %c2(1)
      (-1/Eu1)*P(s_,2)*sigma*c2(2)^(-sigma-1);...                                               %c2(2)
      0; ...                                                                                  %l1(1)
     0; ...                                                                                   %l1(2)
     0; ...                                                                                   %u2btildprime(1)
     0; ...                                                                                   %u2btildprime(2)
     ];
 
 % Derivative of the resource constraint
 dR=[1 0; ...                                                                                    %c1(1)
     0 1; ...                                                                                    %c1(2)
     1 0; ...                                                                                    %c2(1)
     0 1; ...                                                                                    %c2(2)
     -theta_1 0; ...                                                                             %l1(1)    
     0 -theta_1; ...                                                                             %l1(2)
     0 0; ...                                                                                    %u2btildprime(1)    
     0 0;];                                                                                      %u2btildprime(1) 
 
 
grad_eq=[dI dB dR];
end