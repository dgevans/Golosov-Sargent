function [c1,c2,grad_c1,grad_c2] = computeC2_2old(c1_1,c1_2,c2_1,R,s_,P,sigma)
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

    %Compute c2_2 from formula
    frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) ); 
    c2_2 = frac^(-1/sigma); 
    
    
    % vectorize c1, c2 to be of the form described above
    c1=kron(ones(3,1),[c1_1 c1_2]);
    c2=kron(ones(3,1),[c2_1 c2_2]);
    

    
    
    % take gradients
    gradc2_2=zeros(3,1);
    gradc2_2(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2)); % <ok - Anmol>
    gradc2_2(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R; % <ok - Anmol>
    gradc2_2(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2); % <ok - Anmol>
    grad_c1(:,1)=[1;0;0];
    grad_c1(:,2)=[0;1;0];
    grad_c2(:,1)=[0;0;1];
    grad_c2(:,2)=gradc2_2;
    
end
