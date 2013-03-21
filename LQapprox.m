clear all
zbar =  3;
gamma = 2;
sigma = 1;
c1bar = 1.0;
gbar = .3;
bbar = 0;
beta = .9;
alpha1 = .69;
alpha2 = 1-alpha1;

f = @(c2bar,l1bar) [ c2bar - c1bar + l1bar^(gamma)*c1bar^sigma - bbar*(1-beta);...
    c1bar+c2bar+gbar-zbar*l1bar];

x = fsolve(@(x)f(x(1),x(2)),[0.5;0.5]);
c2bar = x(1);
l1bar = x(2);
phi1 = l1bar^(gamma-1)*c1bar^sigma;
rhobar = c2bar^(-sigma)/c1bar^(-sigma);


A = [-1+(sigma*phi1*l1bar)/c1bar, 1,  gamma*phi1, bbar;
    1                         ,  1,  -zbar, 0;
    beta-1                    ,0      , 0     , -c1bar/sigma;
    0                         ,beta-1, 0 , -c2bar/sigma];

b = [alpha1;-alpha1*phi1;alpha2*rhobar;0];

lambdabar = A\b;

omega(1) = -alpha1*sigma/(2*c1bar)-sigma*phi1*l1bar/(2*c1bar)*lambdabar(1)+(1-beta)*(1+sigma)/(2*c1bar)*lambdabar(3);
omega(2) = -alpha2*rhobar*sigma/(2*c2bar) +(1-beta)*(1+sigma)/(2*c2bar)*lambdabar(4);
omega(3) = alpha1*(1-gamma)*phi1/(2*l1bar)-gamma*phi1/(2*l1bar)*lambdabar(1);
omega(4) = lambdabar(1);
omega(5) = lambdabar(1)*phi1*l1bar/2;
omega(6) = -lambdabar(2);
omega(7) = -lambdabar(3);
omega(8) = -lambdabar(4);
omega=omega'
omega(3)+gamma^2/l1bar^2*omega(5)

omegahat(1) = c1bar^2*omega(1);
omegahat(2) = c2bar^2*omega(2);
omegahat(3) = l1bar^2*omega(3);
omegahat(4) = bbar*beta*omega(4);
omegahat(5) = omega(5);
omegahat(6) = zbar*l1bar*omega(6);
omegahat(7) = beta*c1bar*omega(7);
omegahat(8) = beta*c2bar*omega(8);
omegahat = omegahat'