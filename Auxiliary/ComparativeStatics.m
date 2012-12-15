clc
clear all
close all

% This file does the comparative statics of the rest points with technology
% and preference parameters

InitData=load('Data/temp/csigmaLow.mat');
Para=InitData.Para;
theta20=1;
thetaRatioMin=1.5;
thetaRatioMax=5;
thetaRatioGridSize=25;
x0=1;
R0=2;
findSteadyState( x0,R0,Para)
thetaRatioGrid=linspace(thetaRatioMin,thetaRatioMax,thetaRatioGridSize);
for thetaind=1:thetaRatioGridSize
    Para.theta_1=thetaRatioGrid(thetaind)*theta20;
    Para.theta_2=.5;
    S(thetaind)=exp(thetaRatioGrid(thetaind))*10;
[ R(thetaind),x(thetaind),PolicyRules ] = findSteadyState( x0,R0,Para);
x0=x(thetaind);
R0=R(thetaind);
        
            c1=PolicyRules(1:2);
    c2=PolicyRules(3:4);
    l1=PolicyRules(5:6);
    l2=PolicyRules(7:8);
    ul2=(1-Para.psi)./(1-l2);
    uc2=Para.psi./(c2.^(Para.sigma));
    ul1=(1-Para.psi)./(1-l1);
    uc1=Para.psi./(c1.^(Para.sigma));
    
    % TAU - From the WAGE optimality of Agent 2
    Tau(thetaind,:)=1-(ul1./(Para.theta_1.*uc1));
        y(1)=c1(1)*Para.n1+c2(1)*Para.n2+Para.g(1);
    y(2)=c1(2)*Para.n1+c2(2)*Para.n2+Para.g(2);
TaxRevenue(thetaind,:)=y.*Tau(thetaind,:);
    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    TransRate(thetaind,:)=(c2-l2.*ul2./uc2)./y;
    

end

figure()
plot(thetaRatioGrid'*theta20 ,TaxRevenue)
title('tauY')
xlabel('theta_1')
ylabel('Tax REvenues')
figure()
subplot(1,2,1)
plot(thetaRatioGrid'*theta20 ,Tau)
subplot(1,2,2)
plot(thetaRatioGrid'*theta20 ,TransRate)

figure()
scatter(x',R',S','.')
xlabel('x')
ylabel('R')

Para=InitData.Para;
x=[]
R=[]
Para=InitData.Para;
alpha1Min=.15;
alpha1Max=.75;
alpha1GridSize=25;
x0=1;
R0=2;
alpha1Grid=linspace(alpha1Min,alpha1Max,alpha1GridSize);
for alpha1ind=1:alpha1GridSize
    Para.alpha_1=alpha1Grid(alpha1ind);
    Para.alpha_2=1-Para.alpha_1;
    S(alpha1ind)=exp(alpha1Grid(alpha1ind)*10);
[ R(alpha1ind),x(alpha1ind),PolicyRules ] = findSteadyState( x0,R0,Para);
x0=x(alpha1ind);
R0=R(alpha1ind);

            c1=PolicyRules(1:2);
    c2=PolicyRules(3:4);
    l1=PolicyRules(5:6);
    l2=PolicyRules(7:8);
    ul2=(1-Para.psi)./(1-l2);
    uc2=Para.psi./(c2.^(Para.sigma));
    ul1=(1-Para.psi)./(1-l1);
    uc1=Para.psi./(c1.^(Para.sigma));
    
    % TAU - From the WAGE optimality of Agent 2
    Tau(alpha1ind,:)=1-(ul1./(Para.theta_1.*uc1));
        y(1)=c1(1)*Para.n1+c2(1)*Para.n2+Para.g(1);
    y(2)=c1(2)*Para.n1+c2(2)*Para.n2+Para.g(2);
TaxRevenue(alpha1ind,:)=y.*Tau(thetaind,:);

    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    TransRate(alpha1ind,:)=(c2-l2.*ul2./uc2)./y;
    


end


figure()
plot(alpha1Grid,x)
figure()
plot(alpha1Grid,R)

figure()
scatter(x',R',S','.')
xlabel('x')
ylabel('R')