% Set Params


if strcmp(computer,'PCWIN')
    sl='\';
    coresize=2;
else
    sl='/';
    coresize=8;
end


    
% 1. Paramters describing the preferences
theta_1=2; % type of Agent 1
theta_2=1; % Type of Agent 2
psi=.6; % Leisure consumption substitution
beta=.96 ;% subjective time discount factor;
n1=.5;
n2=.5;
% 2. Technology
g_l=.15; % Government expenditure in low state s_l
g_h=.17; % Government expenditure in high state s_h
P=[.5 .5;.5 .5]; % Transition Matrix for g shocks
%P=[.75 .25;.75 .25]; % Transition Matrix for g shocks
alpha_1=.5;
alpha_2=1-alpha_1;
alpha_1=alpha_1*n1;
alpha_2=alpha_2*n2;
alpha=[alpha_1 alpha_2]; % Pareto Weights for agent 1 and Agent 2;
sSize=2; % Dimension of the markov state
% 3. Others
pertub=0.00;
ctol=1e-7;
grelax=.95;
Niter=500;
ResolveCtr=1;
NumSim=10000;
btild_1=0;

  
 ApproxMethod='cheb';
  u2btildGridSize=10;
  RGridSize=10;
  OrderOfAppx_u2btild=5;
  OrderOfApprx_R=5;

   ApproxMethod='spli';
  u2btildGridSize=20;
  RGridSize=20;
  OrderOfAppx_u2btild=10;
  OrderOfApprx_R=10;
 
  pwdd=pwd;
compeconpath=[pwd s1 'compecon2011' sl];
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
Para.u2btildGridSize=u2btildGridSize;
Para.RGridSize=RGridSize;
Para.texpath=texpath;
Para.plotpath=plotpath;
Para.datapath=datapath;
Para.ApproxMethod=ApproxMethod;
Para.OrderOfApprx_R=OrderOfApprx_R;
Para.OrderOfAppx_u2btild= OrderOfAppx_u2btild;
Para.grelax=grelax;
Para.ResolveCtr=ResolveCtr;
Para.NumSim=10000;
Para.btild_1=btild_1;

% Document the table for caliberations
rowLabels = {'$\psi$','$\beta$', '$g_{l}$','$g_{h}$'};
columnLabels = {};
matrix2latex([psi beta g_l g_h]', [texpath 'Param.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
rowLabels = {'$\theta$','$\alpha$','$n$'};
columnLabels = {'Agent 1','Agent 2'};
matrix2latex([theta_1 alpha_1 n1;theta_2 alpha_2 n2]', [texpath 'AgentParam.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
rowLabels={'$g_l$','$g_h$'};
columnLabels={'$g_l$','$g_h$'};
matrix2latex(P, [texpath 'Pr.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');
