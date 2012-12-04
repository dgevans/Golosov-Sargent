function [ res] = ResNecessaryConditionForSS(c2xR,s_,Para) 
                c2=c2xR(1:2);
                u2btild=c2xR(3);
                R=c2xR(4);
                Rprime=[R R];
                  psi = Para.psi;
    sigma = Para.sigma;
    beta =  Para.beta;
    P = Para.P;
    theta_1 = Para.theta_1;
    theta_2 = Para.theta_2;
    g = Para.g;
    alpha_2 = Para.alpha_2;
    n1 = Para.n1;
    n2 = Para.n2;
    c1=R^(1/sigma)*c2;
l2 = ( n1*c1+n2*c2+g+n1*theta_2*Rprime-n1*theta_1  )./(theta_2*(n2+Rprime*n1));
l1 = 1 - (1-l2).*Rprime*theta_2/theta_1;
uc2=psi*c2.^(-sigma);
ucc2=-sigma*psi*c2.^(-sigma-1);
 Euc2 = psi*c2.^(-sigma)*(P(s_,:)');   
Term1=alpha_2*c2.^(-sigma);
Term2=(1-R^(1/sigma))*psi*(1-sigma)*c2.^(-sigma);
Term3=((1-psi)*R./(1-l1).^2)*(1/(theta_2*(n2+R*n1)))*(R*theta_2/theta_1);
Term4=((1-psi)./((1-l2).^2)).*(1/(theta_2*(n2+R*n1)));
Term5=u2btild/beta;
Term6=uc2.*P(s_,:).*ucc2*(-1/(Euc2)^2)+ (1/Euc2)*ucc2;
resTemp=Term1./(Term2+(Term3-Term4)-Term5.*Term6);
res=resTemp(1)-resTemp(2);
end

