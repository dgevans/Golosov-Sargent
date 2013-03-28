 clc
 clear all
 close all
 matlabpool open
run('~/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Main/SetPath')
%load ('~/projects/Golosov-Sargent/Data/temp/cproductivity.mat');
load ('~/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Data/temp/cinequality.mat');
g=Para.g;
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;
sigma=Para.sigma;


% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=0;
s_=1;
c10guess=1;
c20guess=.5;
disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,s_,Para,c,V),[ c10guess c20guess],options);

c10 = x(1);
c20 = x(2);
R0=(c10/c20)^(sigma);
TotalResources=(c10*n1+c20*n2+g);
DenL2=theta_2(s_)*R0*n1+theta_2(s_)*n2;
l20=(TotalResources-theta_1(s_)*n1+ theta_2(s_)*n1*R0)/(DenL2);
l10= 1-(1-l20)*theta_2(s_)/theta_1(s_)*R0;
x0=-(c20-c10)*(psi*c20^(-sigma))-((l10/(1-l10))*R0-l20/(1-l20))*(1-psi)+btild_1*psi*c20^(-sigma);
R0=c20^(-sigma)/c10^(-sigma);



%-- Simulate the MODEL -------------------------------------------------
NumSim=500;

K=100;
ex(1).casename='b_{-1}=-1'; 
Para.saveSimPath= ['~/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Data/temp/BootStrapIneq.mat'];

parfor ctrb=1:K
    rHist0 = rand(NumSim,1);
    CoeffFileName='~/Dropbox/2011RA/FiscalPolicy/OrganizedCode/Golosov-Sargent/Data/temp/cinequality.mat';
    SD(ctrb)=RunSimulationsFromT1AltThetaShocks(CoeffFileName,x0,R0,NumSim,Para,rHist0);
end
save(Para.saveSimPath,'SD');
figure()
for i=1:K
AutoCorr(i)=corr(SD(i).TauHist(3:end),SD(i).TauHist(2:end-1))
%CorrOutput(i)=corr(SD(i).TauHist(1:end),SD(i).YHist(1:end))
end