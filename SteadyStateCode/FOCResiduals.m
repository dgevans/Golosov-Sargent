function [ res ] = FOCResiduals(  X,R,x,VR,Vx,Para )
%FOCRESIDUALS Summary of this function goes here
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
U = Para.U;


if min(theta_1) > 0
    c1 = X(1:2);
    c2 = X(3:4);
    l1 = X(5:6);
    l2 = X(7:8);
    Rprime = X(9:10);
    xprime = X(11:12);
    mu = X(13:14);
    lambda = X(15);
    phi = X(16:17);
    xi = X(18:19);
    rho = X(20:21);

    [~,uc1,ul1,ucc1,ull1] = U(c1,l1);
    [~,uc2,ul2,ucc2,ull2] = U(c2,l2);
    
    Euc2 = dot(P,uc2);
    Emu_uc2 = dot(P,mu.*uc2);



    res(1:2) = x*uc2./(beta*Euc2) - uc2.*(c2-c1) - xprime + Rprime.*l1.*ul1 - l2.*ul2;

    res(3) = dot(P,uc1.*(Rprime-R));

    res(4:5) = theta_2.*ul1.*Rprime - theta_1.*ul2;

    res(6:7) = n1*theta_1*l1 + n2*theta_2*l2 - n1*c1 - n2*c2 - g;

    res(8:9) = uc1.*Rprime - uc2;

    res(10:11) = alpha_1*uc1 + mu.*uc2 + lambda*ucc1.*(Rprime-R) - n1*xi + rho.*ucc1.*Rprime;
    
    res(12:13) = alpha_2*uc2 + mu.*( x*ucc2/(beta*Euc2) - uc2 - ucc2.*(c2-c1))...
        - (x*ucc2*Emu_uc2)/(beta*Euc2^2) - n2*xi - ucc2.*rho;


    res(14:15) = alpha_1*ul1 + mu.*Rprime.*( ul1 + l1.*ull1 )...
        + phi.*theta_2.*Rprime.*ull1 + n1*theta_1.*xi;

    res(16:17) = alpha_2*ul2 - mu.*( ul2 + l2.*ull2 ) - phi.*theta_1.*ull2...
        + n2*theta_2.*xi;

    res(18:19) = beta*Vx - mu;

    res(20:21) = beta*VR + mu.*l1.*ul1 + lambda*uc1 + phi.*theta_2.*ul1 + rho.*uc1;
    
else
    
    c1 = X(1:2);
    c2 = X(3:4);
    l1 = X(5:6);
    Rprime = X(7:8);
    xprime = X(9:10);
    mu = X(11:12);
    lambda = X(13);
    xi = X(14:15);
    rho = X(16:17);
    
    
    [~,uc1,ul1,ucc1,ull1] = U(c1,l1);
    [~,uc2,~,ucc2,~] = U(c2,0.5*ones(2,1));

    Euc2 = dot(P,uc2);
    Emu_uc2 = dot(P,mu.*uc2);



    res(1:2) = x*uc2./(beta*Euc2) - uc2.*(c2-c1) - xprime + Rprime.*l1.*ul1;

    res(3) = dot(P,uc1.*(Rprime-R));

    res(4:5) = n1*theta_1*l1 - n1*c1 - n2*c2 - g;

    res(6:7) = uc1.*Rprime - uc2;

    res(8:9) = alpha_1*uc1 + mu.*uc2 + lambda*ucc1.*(Rprime-R) - n1*xi + rho.*ucc1.*Rprime;
    
    res(10:11) = alpha_2*uc2 + mu.*( x*ucc2/(beta*Euc2) - uc2 - ucc2.*(c2-c1))...
        - (x*ucc2*Emu_uc2)/(beta*Euc2^2) - n2*xi - ucc2.*rho;


    res(12:13) = alpha_1*ul1 + mu.*Rprime.*( ul1 + l1.*ull1 )...
        + n1*theta_1.*xi;

    res(14:15) = beta*Vx - mu;

    res(16:17) = beta*VR + mu.*l1.*ul1 + lambda*uc1 + rho.*uc1;
end

end

