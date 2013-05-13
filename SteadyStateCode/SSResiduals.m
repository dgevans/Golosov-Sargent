function [ res ] = SSResiduals( X,Para )
%SSRESIDUALS Summary of this function goes here
%   Detailed explanation goes here

beta = Para.beta;
alpha_1 = Para.alpha_1;
alpha_2 = Para.alpha_2;
theta_1 = Para.theta_1;
theta_2 = Para.theta_2;
n1=Para.n1;
n2=Para.n2;
g = Para.g;
P = Para.P;
P = P(1,:);
S = length(P);
U = Para.U;

if(min(theta_2) > 0 && min(theta_1) > 0)
    c1 = X(1:S);
    c2 = X(S+1:2*S);
    l1 = X(2*S+1:3*S);
    l2 = X(3*S+1:4*S);
    x = X(4*S+1);
    R = X(4*S+2);
    mu = X(4*S+3);
    lambda = X(4*S+4);
    phi = X(4*S+5:5*S+4);
    xi = X(5*S+5:6*S+4);
    rho = X(6*S+5:7*S+4);

   % [~,uc1,ul1,ucc1,ull1] = U(c1,l1,Para);
    %[~,uc2,ul2,ucc2,ull2] = U(c2,l2,Para);

   [~,uc1,ul1,ucc1,ull1] = U(c1,l1);
   [~,uc2,ul2,ucc2,ull2] = U(c2,l2);

    Euc2 = dot(P,uc2);
    Euc1 = dot(P,uc1);

    res(1:S) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - x +R*l1.*ul1-l2.*ul2;

    res(S+1:2*S) = theta_2.*R.*ul1-theta_1.*ul2;

    res(2*S+1:3*S) = n1*theta_1.*l1+n2*theta_2.*l2-n1*c1-n2*c2-g;

    res(3*S+1:4*S) = uc1*R-uc2;

    res(4*S+1:5*S) = alpha_1*uc1 + mu*uc2 - n1*xi + rho.*ucc1*R;

    res(5*S+1:6*S) = alpha_2*uc2 - mu*( uc2 + ucc2.*(c2-c1) ) - n2*xi - rho.*ucc2; 

    res(6*S+1:7*S) = alpha_1*ul1 + mu*R*( ul1 + l1.*ull1 ) + phi.*ull1.*theta_2*R + n1*theta_1.*xi;

    res(7*S+1:8*S) = alpha_2*ul2 - mu*( ul2 + l2.*ull2 ) - phi.*ull2.*theta_1 + n2*theta_2.*xi;

    res(8*S+1:9*S) = -lambda*beta*Euc1 + mu*l1.*ul1 + lambda*uc1 + phi.*ul1.*theta_2 + rho.*uc1;
elseif min(theta_2) == 0
    c1 = X(1:S);
    c2 = X(S+1:2*S);
    l1 = X(2*S+1:3*S);
    x = X(3*S+1);
    R = X(3*S+2);
    mu = X(3*S+3);
    lambda = X(3*S+4);
    xi = X(3*S+5:4*S+4);
    rho = X(4*S+5:5*S+4);

%    [~,uc1,ul1,ucc1,ull1] = U(c1,l1,Para);
%    [~,uc2,~,ucc2,~] = U(c2,0.5*ones(1,2),Para);

        [~,uc1,ul1,ucc1,ull1] = U(c1,l1);
    [~,uc2,~,ucc2,~] = U(c2,0.5*ones(1,2));

    Euc2 = dot(P,uc2);
    Euc1 = dot(P,uc1);
    
    res(1:S) = x*uc2./(beta*Euc2)-uc2.*(c2-c1) - x +R*l1.*ul1;

    res(S+1:2*S) = n1*theta_1.*l1 - n1*c1 - n2*c2 - g;

    res(2*S+1:3*S) = uc1*R-uc2;

    res(3*S+1:4*S) = alpha_1*uc1 + mu*uc2 - n1*xi + rho.*ucc1*R;

    res(4*S+1:5*S) = alpha_2*uc2 - mu*( uc2 + ucc2.*(c2-c1) ) - n2*xi - rho.*ucc2; 

    res(5*S+1:6*S) = alpha_1*ul1 + mu*R*( ul1+l1.*ull1 ) + n1*theta_1.*xi;

    res(6*S+1:7*S) = -lambda*beta*Euc1 + mu*l1.*ul1 + lambda*uc1  + rho.*uc1;
else
    c1 = X(1:S);
    c2 = X(S+1:2*S);
    l2 = X(2*S+1:3*S);
    x = X(3*S+1);
    R = X(3*S+2);
    mu = X(3*S+3);
    lambda = X(3*S+4);
    xi = X(3*S+5:4*S+4);
    rho = X(4*S+5:5*S+4);

%    [~,uc1,ul1,ucc1,ull1] = U(c1,l1,Para);
%    [~,uc2,~,ucc2,~] = U(c2,0.5*ones(1,2),Para);

    [~,uc1,~,ucc1,~] = U(c1,0.5*ones(1,S));
    [~,uc2,ul2,ucc2,ull2] = U(c2,l2);

    Euc2 = dot(P,uc2);
    Euc1 = dot(P,uc1);
    
    res(1:S) = x*uc2./(Euc2)-uc2.*(c2-c1) - beta.*x -l2.*ul2;

    res(S+1:2*S) = n2*theta_2.*l2 - n1*c1 - n2*c2 - g;

    res(2*S+1:3*S) = uc1*R-uc2;

    res(3*S+1:4*S) = alpha_1*uc1 + mu*uc2 - n1*xi + rho.*ucc1*R;

    res(4*S+1:5*S) = alpha_2*uc2 - mu*( uc2 + ucc2.*(c2-c1) ) - n2*xi - rho.*ucc2; 

    res(5*S+1:6*S) = alpha_2*ul2 - mu*( ul2+l2.*ull2 ) + n2*theta_2.*xi;

    res(6*S+1:7*S) = -lambda*beta.*Euc1  + lambda*uc1  + rho.*uc1;
        
    
    

end

