function [l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2)

    %Compute l1 form formula
    l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1; % < ok - Anmol>
    l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);  % < ok - Anmol>
    l1(1) = l1_1num/l1_1den;  % < ok - Anmol>
    l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2; % <ok - Anmol>
    l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2); % <ok - Anmol>
    l1(2) = l1_2num/l1_2den; % <ok - Anmol>
    
    %compute gradients of l1(1) for c1_1,c1_2,c2_1
    l1grad(1,1) = (l1_1den*(n1-n2*theta_1*c2_1/c1_1^2)+l1_1num*n2*c2_1*theta_1/c1_1^2)/l1_1den^2; %<ok - Anmol>
    l1grad(2,1) = 0;  % <ok - Anmol>
    l1grad(3,1) = (l1_1den*(n2+n2*theta_1/c1_1)-l1_1num*n2*theta_1/c1_1)/l1_1den^2;  % <ok - Anmol>
    
    %compute gradients of l1(1) for c1_1,c1_2,c2_1
    l1grad(1,2) = 0; % <ok - Anmol>
    l1grad(2,2) = (l1_2den*(n1-n2*theta_1*c2_2/c1_2^2)+l1_2num*n2*c2_2*theta_1/c1_2^2)/l1_2den^2; % <ok - Anmol>
    l1grad(3,2) = 0; % <ok - Anmol>
    %use chain rule for c2_2
    d_c2_2 = (l1_2den*(n2+n2*theta_1/c1_2)-l1_2num*n2*theta_1/c1_2)/l1_2den^2; % <ok - Anmol>
    l1grad(:,2) = l1grad(:,2)+d_c2_2*grad_c2_2; % <ok - Anmol>
    
    %compute l2 from formula
    l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1; % <ok - Anmol>
    l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
    l2(1) = l2_1num/l2_1den; % <ok - Anmol>
    l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2; % <ok - Anmol>
    l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2; % <ok - Anmol>
    l2(2) = l2_2num/l2_2den; % <ok - Anmol>
    
    %compute gradients of l2(1) for c1_1,c1_2,c2_1
    l2grad(1,1) = (l2_1den*(n1+n1*theta_2/c2_1)-l2_1num*n1*theta_2/c2_1)/l2_1den^2;  % <ok - Anmol>
    l2grad(2,1) = 0; % <ok - Anmol>
    l2grad(3,1) = (l2_1den*(n2-n1*c1_1*theta_2/c2_1^2)+l2_1num*n1*c1_1*theta_2/c2_1^2)/l2_1den^2; % <ok - Anmol>
    
    %compute gradients of l2(2) for c1_1,c1_2,c2_1
    l2grad(1,2) = 0; % <ok - Anmol>
    l2grad(2,2) = (l2_2den*(n1+n1*theta_2/c2_2)-l2_2num*n1*theta_2/c2_2)/l2_2den^2; % <ok - Anmol>
    l2grad(3,2) = 0; % <ok - Anmol>
    %use chain rule to get the effect of c2_2
    d_c2_2 = (l2_2den*(n2-n1*c1_2*theta_2/c2_2^2)+l2_2num*n1*c1_2*theta_2/c2_2^2)/l2_2den^2; % <ok - Anmol>
    l2grad(:,2) = l2grad(:,2)+d_c2_2*grad_c2_2;
    if theta_2==0
        l2(1)=0;
        l2(2)=0;
        l2grad=zeros(3,2);
    end
end


