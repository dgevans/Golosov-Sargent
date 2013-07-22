 clc
 clear all
 close all
 SetPath

  
theta_1_bar=exp(1.4);
theta_2_bar=1;
e1=1.2/100;
e2=3/100;
theta_1=[theta_1_bar*(1-e1) theta_1_bar*(1+e1)];
theta_2=[theta_2_bar*(1-e2) theta_2_bar*(1+e2)];
alpha_1=0.69;
alpha_2=1-alpha_1;
beta_bar=.95;
rr=7/11;
bb=16/19;
P=[rr 1-rr; 1-bb bb];

g_l_y=.11; % g low
g_h_y=.13; % g high
n1=1;  
n2=1;
tau=.2;
g_Y=mean([g_l_y g_h_y]);
AvfFETarget=.5;
z=fsolve(@(z) GetCalibrationFrischElasticity (z,AvfFETarget,theta_1_bar,theta_2_bar,tau,g_Y,n1,n2), [1 1 ]);
gamma=z(1);
Y=z(2);


% BASELINE GOVERNMENT EXPENDITURE LEVELS
g=g_Y*Y;

% BASELINE PSI
psi=1/(1+gamma);

L2=(alpha_2*g + alpha_1*theta_2 - alpha_2*theta_1 - alpha_2*g*psi + alpha_2*psi*theta_1 + alpha_2*psi*theta_2)/(alpha_1*theta_2 + alpha_2*theta_2);
L1=1-(theta_2/theta_1)*(alpha_1/alpha_2)*(1-L2);
C2=(psi/(1-psi))*theta_2*(1-L2);
C1=(alpha_1/alpha_2)*C2;

for s=1:2
    EUc=sum(P(s,:).*C1.^(-1));
    uc=C1(s)^-1;
    Int(s)=uc/(beta_bar*EUc);
d(s)=(Int(s)*beta_bar-1)*100;
end

    
% 1. Paramters describing the preferences
beta=beta_bar*(1+d);% subjective time discount factor;
Para.sigma=1;
n1=1;
n2=1;
% 2. Technology
%P=[.75 .25;.75 .25]; % Transition Matrix for g shocks
alpha_1=.69;
alpha_2=1-alpha_1;
alpha_1=alpha_1*n1;
alpha_2=alpha_2*n2;
alpha=[alpha_1 alpha_2]; % Pareto Weights for agent 1 and Agent 2;
sSize=2; % Dimension of the markov state
% 3. Others
pertub=0.00;	
ctol=1e-7;
grelax=.9;
Niter=200;
ResolveCtr=1;
NumSim=10000;
btild_1=0;

  
 ApproxMethod='cheb';
  xGridSize=10;
  RGridSize=10;
  OrderOfAppx_x=5;
  OrderOfApprx_R=5;

   ApproxMethod='spli';
  xGridSize=25;
  RGridSize=25;
  OrderOfAppx_x=22;
  OrderOfApprx_R=22;
 
  pwdd=pwd;
compeconpath=[pwd sl 'compecon2011' sl];
knitropath=[pwd sl 'knitro' sl];
texpath= [pwd sl 'Tex' sl] ;

plotpath= [pwd sl 'Graphs' sl] ;

datapath='~/Golosov-Sargent/Data/temp/NewCalibrations/';

mkdir(texpath)
mkdir(plotpath)
mkdir(datapath)
addpath(genpath(compeconpath))
addpath(genpath(knitropath))

Para.ctol=ctol;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.psi=psi;
Para.beta=beta ;
Para.g=g;
Para.P=P;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.n1=n1;
Para.n2=n2;
Para.Niter=Niter;
Para.sSize=sSize;
Para.xGridSize=xGridSize;
Para.RGridSize=RGridSize;
Para.texpath=texpath;
Para.plotpath=plotpath;
Para.datapath=datapath;
Para.ApproxMethod=ApproxMethod;
Para.OrderOfApprx_R=OrderOfApprx_R;
Para.OrderOfAppx_x= OrderOfAppx_x;
Para.grelax=grelax;
Para.ResolveCtr=ResolveCtr;
Para.btild_1=btild_1;
Para.order=3;
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
Para.Niter=200; % MAXIMUM NUMBER OF ITERATION


% flagSetRGrid,flagSetxGrid  TAKES TWO VALUES : 0 IF DEFAULT GRID OR 1 FOR USERDEFINED
% GRID

Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-5;
Para.xMax=5;

% EXPERIMENT 1 : TFPIneq
casename='TFPIneqBetaShocks1';
tempbeta=Para.beta;
Para.beta=mean(Para.beta)*ones(1,length(Para.beta));
Para.beta=[.96 .94]

Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 1;
Para.U=@(c,l) UMix(c,l,Para);
[xSS,RSS]=findSteadyState(0,3,Para);
Para.RMin=RSS*.85;
Para.RMax=RSS*1.1;
MainBellman(Para) 


% EXPERIMENT 2 : TFPIneq
Para.xMin=-5;
Para.xMax=5;

casename='TFPIneqBetaShocks2';
tempbeta=Para.beta;
Para.beta=mean(Para.beta)*ones(1,length(Para.beta));
Para.beta=[.97 .93]

Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 1;
Para.U=@(c,l) UMix(c,l,Para);
try
[xSS,RSS]=findSteadyState(0,3,Para);
catch exception
    tempbeta=Para.beta;
Para.beta=mean(Para.beta)*ones(1,Para.sSize)
Para.U=@(c,l) UMix(c,l,Para);
[xSS,RSS,PolicyRules]=findSteadyState(0,3,Para);
Para.beta=tempbeta;
Para.U=@(c,l) UMix(c,l,Para);
end


Para.RMin=RSS*.9;
Para.RMax=RSS*1.1;
MainBellman(Para) 



% EXPERIMENT 3 : TFPIneq
Para.xMin=-5;
Para.xMax=5;

casename='TFPIneqBetaShocks3';
Para.beta=[.98 .9]

Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 1;
tempbeta=Para.beta;
Para.beta=mean(Para.beta)*ones(1,Para.sSize)
Para.U=@(c,l) UMix(c,l,Para);
[xSS,RSS,PolicyRules]=findSteadyState(0,3,Para);
Para.beta=tempbeta;
Para.U=@(c,l) UMix(c,l,Para);
[xSS,RSS,PolicyRules]=findSteadyState(0,3,Para,PolicyRules);

Para.RMin=RSS*.85;
Para.RMax=RSS*1.1;
MainBellman(Para) 

