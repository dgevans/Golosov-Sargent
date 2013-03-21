clear all
casename='sigmaMed';
load(['Data/temp/' casename '.mat'])
clc
% Set up the para structure  
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
x=fsolve(@(x) GetCalibrationFrischElasticity (x,AvfFETarget,theta_1,theta_2,tau,g_Y,n1,n2), [1 1 ]);
gamma=x(1);
Y=x(2);
g=g_Y*Y;
psi=1/(1+gamma);
beta=.9;
alpha_1=0.69;
alpha_2=1-alpha_1;
Para.n1=n1;
Para.n2=n2;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;
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
Para.datapath=['Data/temp/'];
mkdir(Para.datapath)
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName];
options = optimset('Display','off','TolX',1e-10);
                   
 ApproxMethod='cheb';
  u2btildGridSize=15;
  RGridSize=15;
  OrderOfAppx_u2btild=2;
  OrderOfApprx_R=2;
Para.ApproxMethod=ApproxMethod;
Para.OrderOfApprx_R=OrderOfApprx_R;
Para.OrderOfAppx_u2btild= OrderOfAppx_u2btild;

 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[0 4],options);
u2btild_=xState(1);
R_=xState(2);

casenameQuad=[casename 'Quad'];
Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
Para.Niter=150;
Para.sigma = 2;
RGrid.RMin=R_*.9;
RGrid.RMax=R_*1.1;
NewPh=.5;
Para.P=[1-NewPh NewPh;1-NewPh NewPh];
MainBellman(Para,RGrid) 


ex(1).casename=casenameQuad; 


i=1
Para.datapath='Data/temp/'
 Para.plotpath=['Graphs/' ex(i).casename '/'];
 Para.StoreFileName=['c' ex(i).casename '.mat'];
 Para.flagPlot2PeriodDrifts=0

 GetPlotsForFinalSolution3D(Para)

 GetPlotsForFinalSolution(Para)
 
 GetPlotsForFinalSolution3D(Para)
