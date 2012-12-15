% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetPlotsForFinalSolution3D(BellmanData,SimulationData)
close all
plotpath= BellmanData.Para.plotpath;
Para=BellmanData.Para;
c=BellmanData.c;
V=BellmanData.V;
x_state=BellmanData.x_state;
PolicyRulesStore=BellmanData.PolicyRulesStore;
g=Para.g;
Para.P
sigma=Para.sigma;
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

u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];
if nargin==2
ucbtild_bounds(1)=min(SimulationData.xHist)*(1 -sign(min(SimulationData.xHist))*.2);
ucbtild_bounds(2)=max(SimulationData.xHist)*1.2;
Rbounds(1)=min(SimulationData.RHist)*.8;
Rbounds(2)=max(SimulationData.RHist)*1.2;
end

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

surf(xPlotGrid,RPlotGrid,Values','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
xlabel('x')
ylabel('\rho')
zlabel('V')
title('Value Function')
print(gcf,'-dpng',[plotpath 'ValueFunction3D.png'])



%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$
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
        
            c1=PolicyRules(1:2);
    c2=PolicyRules(3:4);
    l1=PolicyRules(5:6);
    l2=PolicyRules(7:8);
    ul2=(1-psi)./(1-l2);
    uc2=psi./(c2.^(sigma));
    ul1=(1-psi)./(1-l1);
    uc1=psi./(c1.^(sigma));
    Rprime=PolicyRules(end-3:end-2);
    % x' - u_c_2* btildprime
    u2btildprime=PolicyRules(end-1:end);
    % btildprime - x'/u_c2
    btildprime=PolicyRules(9:10);
    
    % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul1./(theta_1.*uc1));
    
    % OUTPUT
    y(1)=c1(1)*n1+c2(1)*n2+g(1);
    y(2)=c1(2)*n1+c2(2)*n2+g(2);
    
    % TRANSFERS
    % These are transfers computed on the assumption that Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2-l2.*ul2./uc2;
    

        
        TauMat(xind,Rind,:)=Tau;
        TransMat(xind,Rind,:)=Trans;
        FOCResMat(xind,Rind)=max(abs(fvec));
        u2BtildePrimeMat(xind,Rind,:)=PolicyRules(end-1:end);
        Deltau2BtildePrimeMat(xind,Rind,:)=PolicyRules(end-1:end)-u2btild;
        BtildePrimeMat(xind,Rind,:)=PolicyRules(end-5:end-4);
        RprimeMat(xind,Rind,:)=PolicyRules(end-3:end-2);
        DeltaRPrimeMat(xind,Rind,:)=PolicyRules(end-3:end-2)-R;
        EDeltaXMat(xind,Rind)=sum(Para.P(s_,:).*squeeze(u2BtildePrimeMat(xind,Rind,:))')-u2btild;
        VDeltaXMat(xind,Rind)=(sum(Para.P(s_,:).*(squeeze((u2BtildePrimeMat(xind,Rind,:)))'-[u2btild u2btild]).^2)-EDeltaXMat(xind,Rind).^2)^0.5;
        EDeltaRMat(xind,Rind)=sum(Para.P(s_,:).*squeeze(RprimeMat(xind,Rind,:))')-R;
        VDeltaRMat(xind,Rind)=(sum(Para.P(s_,:).*(squeeze( (RprimeMat(xind,Rind,:)))'-[R R]).^2)-EDeltaRMat(xind,Rind).^2)^.5;
        xTarget=[u2btild R s_];
    VDirect=funeval(c(s_,:)',V(s_),xTarget(1:2));
    CheckFit(xind,Rind)=(VDirect-V_new)/V_new;
    
        
    end
    



end

figure()
surf(xPlotGrid,RPlotGrid,squeeze(Deltau2BtildePrimeMat(:,:,1))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
hold on
surf(xPlotGrid,RPlotGrid,squeeze(Deltau2BtildePrimeMat(:,:,2))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

xlabel('x')
ylabel('\rho')
print(gcf,'-dpng',[plotpath 'DeltaX3D.png'])


figure()
surf(xPlotGrid,RPlotGrid,squeeze(DeltaRPrimeMat(:,:,1))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
hold on
surf(xPlotGrid,RPlotGrid,squeeze(DeltaRPrimeMat(:,:,2))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

xlabel('x')
ylabel('\rho')
zlabel('\rho(s)-\rho')
print(gcf,'-dpng',[plotpath 'DeltaR3D.png'])



figure()
surf(xPlotGrid,RPlotGrid,FOCResMat')
xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'FOCErrors3D.png'])

figure()
surf(xPlotGrid,RPlotGrid,VDeltaXMat','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
camlight left

xlabel('x')
ylabel('\rho')
zlabel('\sigma(x(s)-x)')
print(gcf,'-dpng',[plotpath 'VDeltaX3D.png'])

figure()
surf(xPlotGrid,RPlotGrid,VDeltaRMat','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

xlabel('x')
ylabel('\rho')
zlabel('\sigma(\rho(s)-\rho)')
print(gcf,'-dpng',[plotpath 'VDeltaR3D.png'])


figure()
surf(xPlotGrid,RPlotGrid,CheckFit','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

xlabel('x')
ylabel('\rho')
zlabel('FitError')
print(gcf,'-dpng',[plotpath 'ChebError3D.png'])


figure()
plot(BellmanData.ErrorInSupNorm)
xlabel('Iteration');
ylabel('Sup Norm Difference');
print(gcf,'-dpng',[plotpath 'CoeffConvergenceSupNorm.png'])
