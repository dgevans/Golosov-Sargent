
clear all
clc

Para.theta_1=1.1;
Para.theta_2=0;
Para.beta=.9;
Para.g=[.1 .2];
Para.P=[.5 .5];
options = optimset('Display','off','TolX',1e-10);
[xR,~,exitflag]=fsolve(@(xR) ResFOCWithXprimeEndowmentEconomy(xR,Para) ,[1 3])
R=xR(2);
x=xR(1);
Check1=R-x-1
[res,c2_ss,tau_ss,Der_tau_ss_x]=ResFOCWithXprimeEndowmentEconomy(xR,Para)
xR=fsolve(@(xR) ResFOCWithXprimeEndowmentEconomy(xR,Para) ,xR)
