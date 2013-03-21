clear all
z1bar = 3.3;
z2bar = 1;
gamma = 1;
sigma = 1;
rhobar = 3;
gbar = .3;
bbar = 1;
beta = .9;
alpha1 = .69;
alpha2 = 1-alpha1;

f = @(c1bar,c2bar,l1bar,l2bar) [ c2bar - c1bar + l1bar^(1+gamma)*c1bar^sigma - l2bar^(1+gamma)*c2bar^sigma - bbar*(1-beta);...
    c1bar+c2bar+gbar-z1bar*l1bar-z2bar*l2bar;
    l1bar^gamma*c1bar^sigma/z1bar - l2bar^gamma*c2bar^sigma/z2bar
    c2bar^(-sigma)/c1bar^(-sigma)-rhobar];

x = fsolve(@(x)f(x(1),x(2),x(3),x(4)),[0.5;0.5;0.5;0.5]);
c1bar = x(1);
c2bar = x(2);
l1bar = x(3);
l2bar = x(4);
phi1 = l1bar^(gamma)*c1bar^sigma;
phi2 = l2bar^gamma*c2bar^sigma;
rhobar = c2bar^(-sigma)/c1bar^(-sigma);
I1 = phi1*l1bar; %FIX THIS
I2 = phi2*l2bar; %FIX This


A = [sigma*I1-c1bar   sigma     c1bar      (1-beta)*sigma        0;
    (1+gamma)*I1      gamma   -z1bar*l1bar         0             0;
    c2bar-sigma*I2   -sigma     c2bar              0      (1-beta)*sigma;
    -(1+gamma)*I2    -gamma   -z2bar*l2bar         0             0;
    bbar*beta            0         0             beta           beta]

b = -[alpha1*c1bar;
     -alpha1*phi1*l1bar;
     alpha2*rhobar*c2bar;
     -alpha2*rhobar*phi2*l2bar;
     0];
 

 lambda = A\b;
 
 omega1(1) = -alpha1*sigma*c1bar + lambda(1)*sigma*(sigma-1)*I1 ...
     + lambda(2)*sigma*(sigma-1) - lambda(4)*(1-beta)*sigma*(1+sigma);
 omega1(2) = -(lambda(1)*sigma*(1+gamma)*I1 + lambda(2)*gamma*sigma);
 omega1(3) = lambda(2)*sigma;
 omega1(4) = -alpha1*gamma*phi1*l1bar + lambda(1)*gamma*(1+gamma)*I1...
     + lambda(2)*gamma*(gamma-1);
 omega1(5) = lambda(2)*gamma + lambda(3)*z1bar*l1bar;
 
 omega2(1) = -alpha2*rhobar*c2bar*sigma - lambda(1)*sigma*(sigma-1)*I2...
     -lambda(2)*sigma*(sigma-1) - lambda(5)*(1-beta)*sigma*(sigma+1);
 omega2(2) = lambda(1)*sigma*(1+gamma)*I2 + lambda(2)*gamma*sigma;
 omega2(3) = -lambda(2)*sigma;
 omega2(4) = -alpha2*rhobar*gamma*phi2*l2bar - lambda(1)*gamma*(1+gamma)*I2...
     - lambda(2)*gamma*(gamma-1);
 omega2(5) = -lambda(2)*gamma + lambda(3)*z2bar*l2bar;
 
 
 omegac1 = omega1(1);
 omegac2 = omega2(1);
 
 theta1 = omega1(2)/omega1(1);
 theta2 = omega2(2)/omega2(1);
 
 eta1 = omega1(3)/omega1(1);
 eta2 = omega2(3)/omega2(1);
 
 omegal1 = omega1(4)-omegac1*theta1^2;
 omegal2 = omega2(4)-omegac2*theta2^2;
 
 omegaz1 = (omega1(5)-omegac1*eta1)/omegal1;
 omegaz2 = (omega2(5)-omegac2*eta2)/omegal2;
 
% A = [sigma         gammma       -sigma    -gamma        0           0;
 %     c1bar      -z1bar*l1bar 
 
 
