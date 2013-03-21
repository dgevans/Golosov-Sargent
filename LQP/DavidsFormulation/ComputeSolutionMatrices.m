clc
LQaproxBig
% Let U=[c1 c2 l1 l2 psi xi]', xi=[g z1 z2]' , lambda=[lambda1 lambda2] 

% FOC : 
% A0 U = A1 mu +A2 xi + A3 lambda- lambda_ /beta

A0=[ sigma              -sigma              gamma                       -gamma                          0           0;
     c1bar              c2bar               -z1bar*l1bar               -z2bar*l2bar                     0           0;
      omegac1           0                   -omegac1*theta1            0                               -sigma      -c1bar;
     0                  omegac2               0                          -omegac2*theta2                sigma      -c2bar;
      -omegac1*theta1    0                   theta1^2*omegac2+omegal1   0                               -gamma      z1bar*l1bar;
    0                   omegac2*theta2     0                           theta1^2*omegac2 + omegal2     gamma        z2bar*l2bar;
    ];

A1=[0;
    0;
    sigma*phi1*l1bar-c1bar;
    c2bar-sigma*phi2*l2bar
    phi1*l1bar*(1+gamma);
    -phi2*l2bar*(1+gamma);
    ];


A2=[0       -1                                      1;
    -gbar   z1bar*l1bar                             z2bar*l2bar;
    0       omegac1*eta1                            0;
    0       0                                       omegac2*eta2;
    0       omegal1*omegaz1-theta1*eta1*omegac1    0;
    0       0                                       omegal2*omegaz2-theta2*eta2*omegac2;];


A3=[0 0;
    0 0;
    1 0;
    0 1;
    0 0;
    0 0];


