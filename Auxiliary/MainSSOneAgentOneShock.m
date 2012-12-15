clear all
theta_1GridSize=5;
theta_1Grid=linspace(2.5,5,theta_1GridSize)
for theta_1GridInd=1:theta_1GridSize
 theta_1=theta_1Grid(theta_1GridInd)
alpha_1=0.5;
u2btildMin=-3;
u2btildMax=0;
RMin=4.5;
RMax=5.5;
tic
[PolicyRules,Para]=SolveSSForOneAgentOneShockEconomy(theta_1,alpha_1,u2btildMin,u2btildMax,RMin,RMax);
 load([Para.datapath '/cVHat.mat'])
 options=optimset('TolX',1e-7,'TolFun',1e-7);
 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[mean(Para.u2bdiffGrid) mean(Para.RGrid)])
 XR(theta_1GridInd,:,:)=xState;
end
figure()
plot(theta_1Grid,squeeze(XR(:,:,1)))
Para.StoreFileName='/cVHat.mat'
Para.flagPlot2PeriodDrifts=0
GetPlotsForFinalSolution(Para)
