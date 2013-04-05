clear all
addpath SIMS

MainLQP



Gamma0 = [I1*sigma-c1bar   0    c2bar-sigma*I1     0     I1*gamma1_I     -I2*gamma2_I     bbar*beta   bbar*beta    0            0            0           0           0           0;
              sigma        0        -sigma         0       gamma1            -(gamma2)        0           0        0            0            0           0           0           0;
              -c1bar       0        -c2bar         0     z1bar*l1bar      z2bar*l2bar         0           0        0            0            0           0           0           0;
              -sigma     sigma        0            0         0                 0              0           1        0            0            0           0           0           0;
                0          0        -sigma        sigma      0                 0              0           1        0            0            0           0           0           0;
           Omegac1c1       0          0            0     Omegac1l1             0              0       Omegac1Q  sigma*I1-c1bar  0          sigma      -c1bar      -sigma         0;
                0          0    Omegac2c2          0         0            Omegac2l2           0       Omegac2Q  c2bar-sigma*I2  0         -sigma      -c2bar         0         -sigma;
           Omegac1l1       0          0            0     Omegal1l1             0              0           0     I1*gamma1_I     0          gamma1     z1bar*l1bar    0           0;
                0          0    Omegac2l2          0         0            Omegal2l2           0           0    -I2*gamma2_I     0         -(gamma2)   z2bar*l2bar    0           0;
                0          0          0            0         0                 0              0       Omegab2Q  beta*bbar  -beta*bbar        0           0           0           0;
           Omegac1Q        0    Omegac2Q           0         0                 0          Omegab2Q        0     bbar*beta       0            0           0           1           1;
                1          0          0            0         0                 0              0           0        0            0            0           0           0           0;
                0          0          1            0         0                 0              0           0        0            0            0           0           0           0;
                0          0          0            0         0                 0              0           0        1            0            0           0           0           0];
 
Gamma0(15,:) =...
     [  sigma*(1-taubar)   0          0            0   gamma1*(1-taubar)       0              0           0        0            0            0           0           0           0];
Gamma0(16,:) =...
     [          0          0   sigma*I2-c2bar      0         0             gamma2_I*I2        0           0        0            0            0           0           0           0];       
Gamma0(17,:) =...
               [0          0          0            0     -z1bar*l1bar     -z2bar*l2bar        0           0        0            0            0           0           0           0];
Gamma0(:,15:17) =zeros(17,3); Gamma0(15,15) = taubar; Gamma0(16,16) = Tbar; Gamma0(17,17) = ybar;
           
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
Psi(15,1) = (1-taubar);
Psi(17,1) = z1bar*l1bar;
Psi(17,2) = z2bar*l2bar;

PI = zeros(17,3);
PI(12:14,:) = eye(3);


[G1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(Gamma0,Gamma1,C,Psi,PI,1/sqrt(beta));
impact
%Now Simulate
y0 = zeros(17,1);
T = 100;
shocks = zeros(3,2);
shocks(3,:) = [-0.01,0.01];
P = 0.5*ones(2);
break;
[ yhatHist shatHist shockhatHist] = SimulateFromT1Linear(y0,G1,impact,shocks,P,T,rHist);
globSol = load('cSmallG.mat');
globSol.Para.saveSimPath = 'SimData.mat';
R0 = rhobar;
x0 = bbar*beta*uc2;
[SD]=RunSimulationsFromT1Alt('cSmallG.mat',x0,R0,T,globSol.Para,rHist);