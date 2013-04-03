clc
clear all
MainLQP
N=zeros(11,11);
N(1,:)=[Omegac1c1 0 Omegac1l1   0 Omegac1Q  0 I1*sigma-c1bar    sigma   c1bar 1/beta 0];
N(2,:)=[0   Omegac2c2   0 Omegac2l2 Omegac2Q    0   c2bar-I2*sigma    -sigma  c2bar  0 1/beta];
N(3,:)=[Omegac1l1 0 Omegal1l1   0   0   0   I1*(1+gamma)    gamma   -z1bar*l1bar 0 0 ];
N(4,:)=[0   Omegac2l2   0   Omegal2l2   0   0   -I2*(1+gamma)   -gamma   -z2bar*l2bar 0 0 ];
N(5,:)=[Omegac1Q Omegac2Q 0 0 0 0 beta*b2bar 0 0 0 0];
N(6,:)=[0 0 0 0 Omegab2Q 0 beta*b2bar 0 0 0 0];
N(7,:)=[I1*sigma-c1bar  -I2*sigma+c2bar     I1*(1+gamma)    -I2*(1+gamma) beta*b2bar beta*b2bar 0   0   0   0   0   ];
N(8,:)=[sigma -sigma gamma -gamma 0 0 0 0 0 0 0];
N(9,:)=[c1bar c2bar -z1bar*l1bar -z2bar*l2bar 0 0 0 0 0 0 0];
N(10,1)=1;
N(10,5)=-1/sigma;
N(11,2)=1;
N(11,5)=-1/sigma;

L=zeros(11,11);
L(1,10)=1;
L(2,11)=1;
L(5,:)=[0   0   0   0   0 -Omegab2Q    0    0   0   -1/sigma    -1/sigma];
L(6,:)=[0   0   0   0   0   0   beta*b2bar 0    0   0   0  ];
L(7,6)=-beta*b2bar;
L(10,:)=[1 0 0 0 0 0 0 0 0 0 0];
L(11,:)=[0 1 0 0 0 0 0 0 0 0 0];

%[U,T]=schur(L\N)
%E=ordeig(T)
%[U,T]=ordschur(U,T,abs(E)<1/sqrt(beta))
%qz(L,N)
[A,B,Q,Z,v]=qz(L,N);
[A,B,Q,Z,v] = qzdiv(1/sqrt(beta),A,B,Q,Z,v)


  %[Lbar,Nbar,W,V] = schur(inv(L)*N)