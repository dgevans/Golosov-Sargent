clc
 clear all
 close all
 SetPath

  

    
% 1. Paramters describing the preferences
theta_1_med = 3.3;
theta_2_med = 1;
ProductivityMultiplier_h=1.02;
ProductivityMultiplier_l=2-ProductivityMultiplier_h;
theta_1=[theta_1_med*ProductivityMultiplier_l theta_1_med*ProductivityMultiplier_h] ; % type of Agent 1
theta_2=[theta_2_med*ProductivityMultiplier_l theta_2_med*ProductivityMultiplier_h] ; % type of Agent 2
psi=.69; % Leisure consumption substitution
beta=.9 ;% subjective time discount factor;
Para.sigma=1;
n1=1;
n2=1;
% 2. Technology
g=.3; % Government expenditure
P=[.5 .5;.5 .5]; % Transition Matrix for g shocks
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
Niter=500;
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
  OrderOfAppx_x=19;
  OrderOfApprx_R=19;
 
  pwdd=pwd;
compeconpath=[pwd sl 'compecon2011' sl];
knitropath=[pwd sl 'knitro' sl];
texpath= [pwd sl 'Tex' sl] ;

plotpath= [pwd sl 'Graphs' sl] ;

datapath=[pwd sl 'Data' sl] ;

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
Para.NumSim=10000;
Para.btild_1=btild_1;
Para.U=@ UMix;

 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
Para.Niter=200; % MAXIMUM NUMBER OF ITERATION


% flagSetRGrid,flagSetxGrid  TAKES TWO VALUES : 0 IF DEFAULT GRID OR 1 FOR USERDEFINED
% GRID

Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-3;
Para.xMax=3;

% EXPERIMENT 1 : productivity
casename='productivity';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.sigma = 1;
Para.RMin=2.7;
Para.RMax=3.2;
MainBellman(Para) 

% --- SOLVE THE BEllMAN FOR INEQUALITY SHOCKS -----------
meantheta = mean([theta_1,theta_2]);
shockSize = 0.0075*2;
%shockSize = 0.0;
Para.theta_1 = [theta_1(1)+shockSize*meantheta theta_1(2)-shockSize*meantheta];
Para.theta_2 = [theta_2(1)-shockSize*meantheta theta_2(2)+shockSize*meantheta];
Para.sigma = 1;
Para.U=@ UMix;
Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-3;
Para.xMax=3;

% EXPERIMENT 2 : inequality
casename='ignequality';
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.RMin=2.7;
Para.RMax=3.1;
%BellmanData=load('~/Data/temp/cTFPIneq.mat')
MainBellman(Para) 


