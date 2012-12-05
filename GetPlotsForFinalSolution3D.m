% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetPlotsForFinalSolution3D(Para,Domain)
flagPlot2PeriodDrifts=Para.flagPlot2PeriodDrifts;
oldplotpath=Para.plotpath;
load([Para.datapath Para.StoreFileName])
close all;
plotpath=oldplotpath;
mkdir(plotpath);
datapath='Data/Calibration/';
disp('Govt Exp')
g=Para.g;
Para.P
n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
%disp('Govt Exp')
%g=Para.g
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;

xSolved=x_state(IndxSolved,:);
xUnSolved=x_state(IndxUnSolved,:);


u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];
if nargin==2
ucbtild_bounds=Domain.xBounds;
Rbounds=Domain.RBounds;
end
%Caption : fig: 3D value functions

figure()
s_=1;
xPlotGridSize=20;
RPlotGridSize=20;
xPlotGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),xPlotGridSize);
RPlotGrid=linspace(Rbounds(1),Rbounds(2),RPlotGridSize);

for xind=1:xPlotGridSize
    for Rind=1:RPlotGridSize
        Values(xind,Rind)=funeval(c(s_,:)',V(s_),[xPlotGrid(xind) RPlotGrid(Rind)]);
    end
end

surf(xPlotGrid,RPlotGrid,Values')
xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'ValueFunction3D.png'])



%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$
s_=1;
s_=1;
xPlotGridSize=20;
RPlotGridSize=20;
xPlotGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),xPlotGridSize);
RPlotGrid=linspace(Rbounds(1),Rbounds(2),RPlotGridSize);

for Rind=1:RPlotGridSize
    for xind=1:xPlotGridSize
        R=RPlotGrid(Rind);
        u2btild=xPlotGrid(xind);
        [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
        if exitflag==1
            IndxPrint(xind)=1;
        else
            IndxPrint(xind)=0;
        end
        
        FOCRes(xind,Rind)=max(abs(fvec));
        u2BtildePrime(xind,Rind,:)=PolicyRules(end-1:end);
        Deltau2BtildePrime(xind,Rind,:)=PolicyRules(end-1:end)-u2btild;
  
        BtildePrime(xind,Rind,:)=PolicyRules(end-5:end-4);
        Rprime(xind,Rind,:)=PolicyRules(end-3:end-2);
        DeltaRPrime(xind,Rind,:)=PolicyRules(end-3:end-2)-R;
        EDeltaX(xind,Rind)=sum(Para.P(s_,:).*squeeze(u2BtildePrime(xind,Rind,:))')-u2btild;
        VDeltaX(xind,Rind)=(sum(Para.P(s_,:).*(squeeze((u2BtildePrime(xind,Rind,:)))'-[u2btild u2btild]).^2)-EDeltaX(xind,Rind).^2)^0.5;
        EDeltaR(xind,Rind)=sum(Para.P(s_,:).*squeeze(Rprime(xind,Rind,:))')-R;
        VDeltaR(xind,Rind)=(sum(Para.P(s_,:).*(squeeze( (Rprime(xind,Rind,:)))'-[R R]).^2)-EDeltaR(xind,Rind).^2)^.5;
        
    end
    



end



figure()
surf(xPlotGrid,RPlotGrid,squeeze(Deltau2BtildePrime(:,:,1)))
hold on
surf(xPlotGrid,RPlotGrid,squeeze(Deltau2BtildePrime(:,:,2)))
xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'DeltaX3D.png'])




figure()
surf(xPlotGrid,RPlotGrid,FOCRes')
xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'FOCErrors3D.png'])

figure()
surf(xPlotGrid,RPlotGrid,VDeltaX','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
camlight left
xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'VDeltaX3D.png'])

figure()
surf(xPlotGrid,RPlotGrid,VDeltaR','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
camlight left
xlabel('x')
ylabel('R')
ylabel('R')
print(gcf,'-dpng',[plotpath 'VDeltaR3D.png'])


