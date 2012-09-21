function [ c2_2 grad ] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma)

    %Compute c2_2 from formula
    frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) ); % <ok - Anmol>
    c2_2 = frac^(-1/sigma); % <ok - Anmol>
    grad=zeros(3,1);
    %compute the gradients for c1_1,c1_2,c2_1
    grad(1) = c1_1^(-sigma-1)*frac^(-1/sigma-1)*R*P(s_,1)/(P(s_,2)); % <ok - Anmol>
    grad(2) = c1_2^(-sigma-1)*frac^(-1/sigma-1)*R; % <ok - Anmol>
    grad(3) = -c2_1^(-sigma-1)*frac^(-1/sigma-1)*P(s_,1)/P(s_,2); % <ok - Anmol>
end

