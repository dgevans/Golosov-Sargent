% this file solves the bellman equation in BEGS for alternative pareto wts 

clc
 clear all
 close all
 SetPath

 SetParaStruc
theta_1=1/3; % theta low
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
alphaMin=0.6;
alphaMax=.67;




alphaGrid=linspace(alphaMin,alphaMax,gridSize);
Para.U = @(c,l) UMix(c,l,Para);
x=0;
R=1/3;
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
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
Para.Niter=200; % MAXIMUM NUMBER OF ITERATION


% case 1

casename='LowAlpha1';
alpha_1=.05;
alpha_2=1-alpha_1;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-1;
Para.xMax=1;
Para.U=@(c,l) UMix(c,l,Para);


Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
[~,indx]=min(( reversealpha1grid-alpha_1).^2);
Para.RMin=RR(indx)*.8;
Para.RMax=RR(indx)*1.2;
[ xSS,RSS,PolicyRule ] = findSteadyState( X(indx),RR(indx),Para,PolicyRules(indx,:));
Para.xSS=xSS;
Para.RSS=RSS;
%MainBellman(Para) 

% case 2


casename='MidAlpha1';
alpha_1=.5;
alpha_2=1-alpha_1;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-1;
Para.xMax=1;
Para.U=@(c,l) UMix(c,l,Para);


Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
[~,indx]=min(( reversealpha1grid-alpha_1).^2);
Para.RMin=RR(indx)*.75;
Para.RMax=RR(indx)*1.25;
[ xSS,RSS,PolicyRule ] = findSteadyState( X(indx),RR(indx),Para,PolicyRules(indx,:));
Para.xSS=xSS;
Para.RSS=RSS;
%MainBellman(Para) 

% case 2

casename='HighAlpha1';
alpha_1=.8;
alpha_2=1-alpha_1;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.flagSetRGrid=1; 
Para.flagSetxGrid=1;
Para.xMin=-.7;
Para.xMax=.7;
Para.U=@(c,l) UMix(c,l,Para);

Para.StoreFileName=['c' casename '.mat'];
CoeffFileName=[Para.datapath Para.StoreFileName]; 
[~,indx]=min(( reversealpha1grid-alpha_1).^2);
Para.RMin=RR(indx)*.8;
Para.RMax=RR(indx)*1.1;
[ xSS,RSS,PolicyRule ] = findSteadyState( X(indx),RR(indx),Para,PolicyRules(indx,:));
Para.xSS=xSS;
Para.RSS=RSS;
%MainBellman(Para) 

