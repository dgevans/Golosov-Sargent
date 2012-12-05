


% Given x,R
% steady state tax rates

function [res,c2_ss,tau_ss,Der_tau_ss_x]= ResFOCWithXprimeEndowmentEconomy(xR,Para)
x=xR(1);
R=xR(2);
theta_1=Para.theta_1;
theta_2=Para.theta_2;
beta=Para.beta;
g=Para.g;
P=Para.P;

TotalConsumption=theta_1+theta_2-g;
c2_ss=TotalConsumption/(1+R);
Ec2_ssinv=sum((1./c2_ss.*P));
Ec_ssinv=sum((1./TotalConsumption.*P));

tau_ss = 1- (R-1+x./(beta*c2_ss.*Ec2_ssinv)-x).*c2_ss/(theta_1-theta_2);

Etau_ss=sum(tau_ss.*P);


Der_tau_ss_x=TotalConsumption./((1+R)*(theta_1-theta_2));


Num=x*beta*Etau_ss;
Den=(1+beta)*beta*(1+R)*(theta_1-theta_2)*Ec_ssinv.*tau_ss;
res=Der_tau_ss_x-Num./Den;

end
