
% This script creates the Para structure that stores the information
% required for solving the value function
    
% 1. Paramters describing the preferences
theta_1=2; % type of Agent 1
theta_2=1; % Type of Agent 2
psi=.6; % Leisure consumption substitution
sigma = 2; % risk aversion
beta=.96 ;% subjective time discount factor;
n1=.5; % mass of type 1 agent
n2=.5; %mass of type 2 agent
% 2. Technology
g_l=.15; % Government expenditure in low state s_l
g_h=.17; % Government expenditure in high state s_h
P=[.5 .5;.5 .5]; % Transition Matrix for g shocks
alpha_1=.5; % pareto wt of agent 1
alpha_2=1-alpha_1; % pareto wt of agent 2
alpha_1=alpha_1*n1; % population adjusted pareto wts
alpha_2=alpha_2*n2; % population adjusted pareto wts
alpha=[alpha_1 alpha_2]; % Pareto Weights for agent 1 and Agent 2;
sSize=2; % Dimension of the markov state
% 3. Others
ctol=1e-8;  % stopping criteria for iner optimization
grelax=.95; % wt of the new coeff
Niter=500; % number of value function iterations
ResolveCtr=1; % Frequency with which the routine for unresolved points must be tried
NumSim=10000; % Number of simulations
btild_1=0; % initial condition for Time0 problem
DeltaX=1; % deviation from steadys state for the default grid
DeltaR=0.5;
%. Polynomial Approximation details
   ApproxMethod='spli'; % basis poly
  xGridSize=15; % density of grid in dimension x
  RGridSize=15; % density of grid in dimension R
  OrderOfAppx_x=10; % number of splines in dimension x
  OrderOfApprx_R=10; % number of splines in dimension R
 
  

% POPULATE THE PARA STRUCTURE
Para.ctol=ctol;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.psi=psi;
Para.sigma = sigma;
Para.beta=beta ;
Para.g_l=g_l;
Para.g_h=g_h;
g=[g_l g_h];
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
Para.DeltaX=DeltaX;
Para.DeltaR=DeltaR;

%DOCUMENT THE PARAMETERS IN A TEX FILE
rowLabels = {'$\psi$','$\beta$', '$g_{l}$','$g_{h}$'};
columnLabels = {};
matrix2latex([psi beta g_l g_h]', [texpath 'Param.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
rowLabels = {'$\theta$','$\alpha$','$n$'};
columnLabels = {'Agent 1','Agent 2'};
matrix2latex([theta_1 alpha_1 n1;theta_2 alpha_2 n2]', [texpath 'AgentParam.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
rowLabels={'$g_l$','$g_h$'};
columnLabels={'$g_l$','$g_h$'};
matrix2latex(P, [texpath 'Pr.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
