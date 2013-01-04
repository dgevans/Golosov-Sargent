function [c1,c2,grad_c1,grad_c2] = computeC2_2(c1,c2_,R,s_,P,sigma)
%COMPUTEC2_@ Computes c_2(2) and return c1 and c2 as 3x2 matrix.  Here c1 will be of
%the form
%     c_1(1)   c_1(2)
%     c_1(1)   c_1(2)
%     c_1(1)   c_1(2)
%Similarly for c2.  The 3x2 matrix format is useful for computing future
%gradients.  gradc1 will a matrix containing the derivative of
%c1 with respect to the various z's.  For instance the ith row and jth
%column will be the derivative of c_1(j) with respect to x(i).  Thus
%gradc1 will be
%       1   0
%       0   1      
%       0   0
%Similarly for gradc2

    S = length(P(1,:));
    P_ = P(s_,:); P_(S) = [];
    %Compute c2_2 from formula
  
    frac = (R*dot(P(s_,:),c1.^(-sigma)) - dot(P_,c2_.^(-sigma)))/P(s_,S);
    
    %frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
    %    /( P(s_,2) ); 
    c2_S = frac^(-1/sigma); 
    
    
   
    
    
    % take gradients
    gradc2_S=zeros(2*S-1,1);
    gradc2_S(1:S,1) = c1.^(-sigma-1)*frac^(-1/sigma-1)*R.*P(s_,:)/P(s_,S);
    %gradc2_2(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2)); % <ok - Anmol>
    %gradc2_2(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R; % <ok - Anmol>
    gradc2_S(S+1:2*S-1,1) = -c2_.^(-sigma-1)*frac^(-1/sigma-1).*P_/P(s_,S);
    
    %gradc2_2(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2); % <ok - Anmol>
    %grad_c1(:,1)=[1;0;0];
    %grad_c1(:,2)=[0;1;0];
    grad_c1 = [eye(S);zeros(S-1,S)];
    %grad_c2(:,1)=[0;0;1];
    %grad_c2(:,2)=gradc2_2;
    grad_c2 = [zeros(S,S-1);eye(S-1)];
    grad_c2(:,S) = gradc2_S;
     % vectorize c1, c2 to be of the form described above
    c1=kron(ones(2*S-1,1),c1);
    c2=kron(ones(2*S-1,1),[c2_ c2_S]);
    

    
end

