% this file solves the bellman equation in BEGS for alternative pareto wts 

clc
 clear all
 close all
 SetPath

 SetParaStruc
theta_1=3; % theta low
theta_2=1;  % theta high
n1=1;  
n2=1;


% BASELINE GOVERNMENT EXPENDITURE LEVELS
g=[.15 .17];

% BASELINE PSI
psi=.69;
% BASELINE DISCOUNT FACTOR

beta=.9;
sigma=1;
% BASELINE PARETO WTS
alpha_1=0.31;
alpha_2=1-alpha_1;
Para.sigma=sigma;
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
Para.g=g;
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
 


gridSize=20;
alphaMin=0.2;
alphaMax=.8;




alphaGrid=linspace(alphaMin,alphaMax,gridSize);
Para.U = @(c,l) UMix(c,l,Para);
x=0;
R=3;
for i=1:length(alphaGrid)
    tic
    Para.alpha_1=alphaGrid(length(alphaGrid)-i+1);
    Para.alpha_2=1-alphaGrid(length(alphaGrid)-i+1);
    Para.U = @(c,l) UMix(c,l,Para);
    if i>1
    [ x,R,PolicyRule ] = findSteadyState( x,R,Para,PolicyRule);
    else
    [ x,R,PolicyRule ] = findSteadyState( x,R,Para);
    end    
    
    PolicyRules(i,:)=PolicyRule;
    X(i)=x;
    RR(i)=R;
    reversealpha1grid(i)=Para.alpha_1;
    toc
end
[reversealpha1grid' X' RR']
