 clc
 clear all
 close all
 SetPath

 %{ 
This file solves the G-S economy with BGP preferences of the form
 psi.c^(1-sigma)/1-sigma+(1-psi)log[1-l] with following calibrations

1.] The ratio of productivities is 3.3 and the low productivity is
 normalized to 1
2.] psi is choosen to get a average FE of about 0.5
3.] pareto wts are such that the no-shock problem gives a tax rate of about
    20 percent
4.] Government expenditures are about 11 and 13 percent of the output
5.] beta =0.9
%}

% - XXXXXXX ANMOL - CURRENTLY IT USES A LEGACY METHOD FOR GETTING THE
% BASELINE PARAMETERS. NOW THAT WE HAVE A STEADY STATE CODE WE CAN USE THIS
% TO TARGET SOME AGREED MOMENTS IN OBSERVABLES 
SetParaStruc
theta_1=3.3; % theta high
theta_2=1;  % theta low
g_l_y=.11; % g low
g_h_y=.13; % g high
n1=1;  
n2=1;
tau=.2;
g_Y=mean([g_l_y g_h_y]);
AvfFETarget=.5;
z=fsolve(@(z) GetCalibrationFrischElasticity (z,AvfFETarget,theta_1,theta_2,tau,g_Y,n1,n2), [1 1 ]);
gamma=z(1);
Y=z(2);


% BASELINE GOVERNMENT EXPENDITURE LEVELS
g=g_Y*Y;

% BASELINE PSI
psi=1/(1+gamma);
% BASELINE DISCOUNT FACTOR

beta=.9;

% BASELINE PARETO WTS
alpha_1=0.69;
alpha_2=1-alpha_1;
Para.n1=n1;
Para.n2=n2;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;

% BASELINE PROBABILITY MATRIX
NewPh=.5;
Para.P=[1-NewPh NewPh;1-NewPh NewPh];

% POPULATE THE PARA STRUC WITH THE BASELINE VALUES
Para.beta=.9;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.psi=psi;
Para.g=[g_l_y g_h_y]*Y;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.btild_1=0;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.datapath=[rootDir sl 'Data/temp/'];
mkdir(Para.datapath)
casename='sigma';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
Para.Niter=200; % MAXIMUM NUMBER OF ITERATION


% flagSetRGrid,flagSetxGrid  TAKES TWO VALUES : 0 IF DEFAULT GRID OR 1 FOR USERDEFINED
% GRID

Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-3;
Para.xMax=3;

% EXPERIMENT 1 : SIGMA=1
casename='sigmaLow';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 1;
Para.RMin=2.2;
Para.RMax=3.5;
MainBellman(Para) 


