clear all

NumSim=5000;
rHist0=rand(NumSim,1);
btild_1=-2;
s_=1;
InitialConditions.s0=s_;



InitialConditions.s0=1;
load ~/Golosov-Sargent/Data/Draft/cTFP.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(1)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneq.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(2)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneqBetaShocks.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(3)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneqLargerBetaShocks.mat
[ InitialConditions.x0,InitialConditions.R0] = solveTime0Problem(Para,c,V,btild_1,s_);
[SimData(4)]=RunSimulationsFromT1AltThetaShocksUsingPolicyRules(PolicyFunctions,InitialConditions,NumSim,Para,rHist0);
save ('~/Golosov-Sargent/Data/Draft/Simulations.mat','SimData','rHist0')
break;


load ~/Golosov-Sargent/Data/Draft/Simulations.mat
a=2;
b=10000;
C={'g','r','b','k'}
%C={'-.k',':k','--k','k'}

figure()
subplot(2,1,1)
for i=1:4
plot(-SimData(i).btildHist(a:end-b)./SimData(i).YHist(a:end-b),C{i},'LineWidth',3)
hold on
end
xlabel('time')
ylabel('Debt-Gdp')




load ~/Golosov-Sargent/Data/Draft/Simulations.mat
a=2;
b=10000;
C={'g','r','b','k'}
%C={'-.k',':k','--k','k'}



%legend('TFP','TFP+Ineq','TFP+Ineq+(small) Beta','TFP+Ineq+(large) Beta')

subplot(2,1,2)
for i=1:4
plot(SimData(i).TauHist(a:end-b),C{i})
hold on
end
axis([1 length(SimData(i).TauHist(a:end-b)) min(SimData(i).TauHist(a:end-b))-3/100 max(SimData(i).TauHist(a:end-b))+3/100])
xlabel('time')
ylabel('Labor Tax rate')



figure()
subplot(2,1,1)
plot([SimData.xHist])
xlabel('time')
ylabel('x_t')
title('Marginal utility adjusted difference in debt')

subplot(2,1,2)
plot([SimData.RHist])
xlabel('time')
ylabel('\rho_t')
title('Ratio of marginal utility')

C={'s','d','+','o'}
fr=30
figure()
%subplot(2,1,1)
for i=1:4
scatter((a:fr:4999)',-SimData(i).btildHist(a:fr:end-b)./SimData(i).YHist(a:fr:end-b),40,'k',C{i},'LineWidth',2)
hold on
end
xlabel('time')
ylabel('Debt-Gdp')

set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')



subplot(2,1,2)
for i=1:4
%plot(SimData(i).TauHist(a:end-b),C{i})
scatter((a:fr:4999)',SimData(i).TauHist(a:fr:end-b),'k',C{i},'LineWidth',1)
hold on
end
axis([1 length(SimData(i).TauHist(a:end-b)) min(SimData(i).TauHist(a:end-b))-3/100 max(SimData(i).TauHist(a:end-b))+3/100])
xlabel('time')
ylabel('Labor Tax rate')

