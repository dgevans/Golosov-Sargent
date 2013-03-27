
L(1,:)=[Omegac1c1 0 Omegac1l1   0 Omegac1Q  0 I1*sigma-c1bar    sigma   c1bar -1 0];
L(2,:)=[0   Omegac2c2   0 Omegac2l2 Omegac2Q    0   c2bar-I2*sigma    -sigma  c2bar  0 -1];
L(3,:)=[Omegac1l1 0 Omegal1l1   0   0   0   I1*(1+gamma)    gamma   -z1bar*l1bar 0 0 ];
L(4,:)=[0   Omegac2l2   0   Omegal2l2   0   0   -I2*(1+gamma)   -gamma   -z2bar*l2bar 0 0 ];
L(5,:)=[Omegac1Q    Omegac2Q    0   0   0 Omegab2Q    beta*b2bar  0   0   1/sigma     1/sigma];
L(6,:)=[0   0   0   0   0   0   beta*b2bar 0    0   0   0  ];
L(7,:)=[I1*sigma-c1bar  -I2*sigma+c2bar     I1*(1+gamma)    -I2*(1+gamma) beta*b2bar beta*b2bar 0   0   0   0   0   ];
L(8,:)=[sigma -sigma gamma -gamma 0 0 0 0 0 0 0];
L(9,:)=[c1bar c2bar -z1bar*l1bar -z2bar*l2bar 0 0 0 0 0 0 0];
L(10,:)=[1 0 0 0 0 0 0 0 0 0 0];
L(11,:)=[0 1 0 0 0 0 0 0 0 0 0];

det(L)

N=zeros(11,11);
N(1,10)=1/beta;
N(2,11)=1/beta;
N(6,5)=Omegab2Q;
N(6,7)=beta*b2bar;
N(7,6)=-beta*b2bar;
N(10,1)=1;
N(10,5)=-1/sigma;
N(11,2)=1;
N(11,5)=-1/sigma;

