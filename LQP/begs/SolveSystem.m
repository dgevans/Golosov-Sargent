clear all
addpath SIMS

MainLQP

Gamma0 = [I1*sigma*c1bar   0    c2bar-sigma*I1     0     I1*(1+gamma)    -I2*(1+gamma)    bbar*beta   bbar*beta    0            0            0           0           0           0;
              sigma        0        -sigma         0       gamma         -(gamma)         0           0        0            0            0           0           0           0;
              -c1bar       0        -c2bar         0     z1bar*l1bar      z2bar*l2bar         0           0        0            0            0           0           0           0;
              -sigma     sigma        0            0         0                 0              0           1        0            0            0           0           0           0;
                0          0        -sigma        sigma      0                 0              0           1        0            0            0           0           0           0;
           Omegac1c1       0          0            0     Omegac1l1             0              0       Omegac1Q  sigma*I1-c1bar  0          sigma      -c1bar      -sigma         0;
                0          0    Omegac2c2          0         0            Omegac2l2           0       Omegac2Q  c2bar-sigma*I2  0         -sigma      -c2bar         0         -sigma;
           Omegac1l1       0          0            0     Omegal1l1             0              0           0     I1*(1+gamma)    0         gamma     z1bar*l1bar    0           0;
                0          0    Omegac2l2          0         0            Omegal2l2           0           0        0       -I2*(1+gamma) -(gamma)   z2bar*l2bar    0           0;
                0          0          0            0         0                 0              0       Omegab2Q  beta*bbar  -beta*bbar        0           0           0           0;
           Omegac1Q        0    Omegac2Q           0         0                 0          Omegab2Q        0        beta*bbar            0            0           0           1           1;
                1          0          0            0         0                 0              0           0        0            0            0           0           0           0;
                0          0          1            0         0                 0              0           0        0            0            0           0           0           0;
                0          0          0            0         0                 0              0           0        1            0            0           0           0           0];
 
Gamma0(15,:) = -(1-taubar)/taubar *...
               [-sigma     0          0            0      -gamma               0              0           0        0            0            0           0           0           0];
Gamma0(16,:) =-[0          0 (c2bar-sigma*I2)/Tbar 0         0           -I2*(1+gamma)/Tbar   0           0        0            0            0           0           0           0];       
Gamma0(17,:) = -1/ybar *...
               [0          0          0            0     z1bar*l1bar       z2bar*l2bar        0           0        0            0            0           0           0           0];
Gamma0(:,15:17) = [zeros(14,3);eye(3)];
           
Gamma1 = zeros(14,14);
Gamma1(1,7) = bbar;
Gamma1(6,13) = -sigma*beta^(-1);
Gamma1(7,14) = -sigma*beta^(-1);
Gamma1(12,2) = 1;
Gamma1(13,4) = 1;
Gamma1(14,10)= 1;
Gamma1(15:17,:) = 0;
Gamma1(:,15:17) = 0;


C = zeros(17,1);

Psi = zeros(17,3);
Psi(2,:) = -[-1,1,0];
Psi(3,:) = -[z1bar*l1bar,z2bar*l2bar,-gbar];
Psi(6,1) = -Omegac1z1;
Psi(7,2) = -Omegac2z2;
Psi(8,1) = -Omegal1z1;
Psi(9,2) = -Omegal2z2;
Psi(15,1) = (1-taubar)/taubar;
Psi(17,1) = z1bar*l1bar/ybar;
Psi(17,2) = z2bar*l2bar/ybar;

PI = zeros(17,3);
PI(12:14,:) = eye(3);


[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(Gamma0,Gamma1,C,Psi,PI);
save('Gamma0','Gamma1','impact','G1')