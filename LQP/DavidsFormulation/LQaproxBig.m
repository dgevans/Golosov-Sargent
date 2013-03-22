clear all
S = 2;
P = 0.5*ones(2,2);
z1bar = 3.3;
z1 = [0;0];
z2bar = 1;
z2 = [0;0];
gamma = 1;
sigma = 1;
rhobar = 2.9;
gbar = .3;
g = [-.05;+.05];
bbar = -3;
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
taubar = 1- l1bar^gamma*c1bar^sigma/z1bar;
Tbar = c2bar-phi2*l2bar;
ybar = z1bar*l1bar+z2bar*l2bar;


A = [sigma*I1-c1bar   sigma     c1bar      (1-beta)*sigma        0;
    (1+gamma)*I1      gamma   -z1bar*l1bar         0             0;
    c2bar-sigma*I2   -sigma     c2bar              0      (1-beta)*sigma;
    -(1+gamma)*I2    -gamma   -z2bar*l2bar         0             0;
    bbar*beta            0         0             beta           beta];

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
 
 thetaz1 = (omega1(5)-omegac1*eta1)/omegal1;
 thetaz2 = (omega2(5)-omegac2*eta2)/omegal2;
 
 A2 = [sigma           gamma                -sigma          -gamma         0           0;
     c1bar         -z1bar*l1bar             c2bar       -z2bar*l2bar       0           0;
     omegac1      -omegac1*theta1             0               0         -sigma       -c1bar;
  -theta1*omegac1 theta1^2*omegac1+omegal1   0               0         -gamma    z1bar*l1bar;
       0                 0                omegac2      -theta2*omegac2   sigma       -c2bar;
       0                 0          -theta2*omegac2  theta2^2*omegac2+omegal2   gamma  z2bar*l2bar;];

%Weightings on each of the multiplier and shock terms
muF = A2\[0;
        0;
        sigma*I1-c1bar;
        (1+gamma)*I1;
        c2bar-sigma*I2;
        -(1+gamma)*I2];
    
z1F = A2\[1;
          z1bar*l1bar;
          omegac1*eta1;
          thetaz1-theta1*eta1*omegac1;
          0;
          0];
      
z2F = A2\[-1;
          z2bar*l1bar;
          0;
          0;
          omegac2*eta2;
          thetaz2-theta2*eta2*omegac2];
gF = A2\[0;
         -gbar;
         0;
         0;
         0;
         0];
     
lambda1F = A2\[0;
               0;
               1;
               0;
               0;
               0];
           
lambda2F = A2\[0;
               0;
               0;
               0;
               1;
               0];
lambdaF = [lambda1F lambda2F];

%Discounted sums of each of the shocks
sumbetaz1 = (eye(S)-beta*P)\z1;
sumbetaz2 = (eye(S)-beta*P)\z2;
sumbetag = (eye(S)-beta*P)\g;
shocks = [z1 z2 g];
sumbetashocks = [sumbetaz1 sumbetaz2 sumbetag];

%Weightings for the taxes transfers and output on the shocks mu and the
%lambdas
mult = [z1F z2F gF];
tauMult = (1-taubar)/taubar *([1 0 0]-gamma*mult(2,:)-sigma*mult(1,:));
muFtau = (1-taubar)/taubar *(-gamma*muF(2)-sigma*muF(1));
lambdaFtau = (1-taubar)/taubar *(-gamma*lambdaF(2,:)-sigma*lambdaF(1,:));
TMult = ((c2bar-I2*sigma)*mult(3,:)-(1+gamma)*I2*mult(4,:))/Tbar;
muFT = ((c2bar-I2*sigma)*muF(3)-(1+gamma)*I2*muF(4))/Tbar;
lambdaFT = ((c2bar-I2*sigma)*lambdaF(3,:)-(1+gamma)*I2*lambdaF(4,:))/Tbar;
yMult = ( z1bar*l1bar*([1,0,0]+mult(2,:))+z2bar*l2bar*([0,1,0]+mult(4,:)))/ybar;
muFy = (z1bar*l1bar*muF(2)+z2bar*l2bar*muF(4))/ybar;
lambdaFy = (z1bar*l1bar*lambdaF(2,:)+z2bar*l2bar*lambdaF(4,:))/ybar;

fMult = (1-beta)*( taubar*ybar*(tauMult+yMult)-2*Tbar*TMult-gbar*[0,0,1]);
muFf = 2*Tbar*muFT-taubar*ybar*(muFtau+muFy);

f = (fMult*sumbetashocks')';
Ef = P*f;
eps = f-Ef;

epsLamb1 = (1-beta)*(mult(1,:)*sumbetashocks')'-(mult(1,:)*shocks')';
epsLamb2 = (1-beta)*(mult(3,:)*sumbetashocks')'-(mult(3,:)*shocks')';

lambMat = [lambda1F(1) lambda2F(1);
           lambda1F(3) lambda2F(3)];
lambWeigts = inv(lambMat); %Weights of the epsLam1 for each of the lambda multipliers
%1st row is for lambda 1, second is for lambda 2



 
 
