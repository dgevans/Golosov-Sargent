% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetTransitionDynamics(BellmanData,SimulationData)
close all
plotpath=BellmanData.Para.plotpath;
Para=BellmanData.Para;
c=BellmanData.c;
V=BellmanData.V;
domain=BellmanData.domain;
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

xLL=Para.xMin;
xUL=Para.xMax;
ucbtild_bounds = [xLL,xUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];
if nargin==2
ucbtild_bounds(1)=min(SimulationData.xHist)*(1 -sign(min(SimulationData.xHist))*.2);
ucbtild_bounds(2)=max(SimulationData.xHist)*1.2;
Rbounds(1)=min(SimulationData.RHist)*.8;
Rbounds(2)=max(SimulationData.RHist)*1.2;
end

s_=1;
xPlotGridSize=20;
RPlotGridSize=20;
xPlotGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),xPlotGridSize);
RPlotGrid=linspace(Rbounds(1),Rbounds(2),RPlotGridSize);

%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$
s_=1;
xPlotGridSize=25;
RPlotGridSize=25;
xPlotGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),xPlotGridSize);
RPlotGrid=linspace(Rbounds(1),Rbounds(2),RPlotGridSize);
ctr=1;
for Rind=1:RPlotGridSize
    for xind=1:xPlotGridSize
        
        R=RPlotGrid(Rind);
        x=xPlotGrid(xind);
                ApproximationDomain(ctr,:)=[xPlotGrid(xind) RPlotGrid(Rind)];

        [PolicyRulesInit]=GetInitialApproxPolicy([x R s_] ,domain,PolicyRulesStore);
        [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(x,R,s_,c,V,PolicyRulesInit,Para);
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
    xprime=PolicyRules(end-1:end);
    % btildprime - x'/u_c2
    btildprime=PolicyRules(9:10);
    
    % TAU - From the WAGE optimality of Agent 2
    Tau=1-(ul1./(theta_1.*uc1));
    
    % OUTPUT
    y(1)=c1(1)*n1+c2(1)*n2+g(1);
    y(2)=c1(2)*n1+c2(2)*n2+g(2);
    
    % TRANSFERS
    % These are transfers computed on the assumption tPlot Agent 2 cannot
    % borrow and lend. The transfers are the difference between his
    % consumption and after tax earning (l . U_l/U_c)
    Trans=c2-l2.*ul2./uc2;
    

        
        TauMat(xind,Rind,:)=sum(Para.P(s_,:).*Tau);
        TransMat(xind,Rind,:)=sum(Para.P(s_,:).*Trans);
        FOCResMat(xind,Rind)=max(abs(fvec));
        xePrimeMat(xind,Rind,:)=PolicyRules(end-1:end);
        xPlotprime(ctr,:)=PolicyRules(end-1:end);
        
        DeltaxePrimeMat(xind,Rind,:)=PolicyRules(end-1:end)-x;
        BtildePrimeMat(xind,Rind,:)=sum(Para.P(s_,:).*PolicyRules(end-5:end-4));
        RprimeMat(xind,Rind,:)=PolicyRules(end-3:end-2);
        RPlotprime(ctr,:)=PolicyRules(end-3:end-2);
        DeltaRPrimeMat(xind,Rind,:)=PolicyRules(end-3:end-2)-R;
        EDeltaXMat(xind,Rind)=sum(Para.P(s_,:).*squeeze(xePrimeMat(xind,Rind,:))')-x;
        VDeltaXMat(xind,Rind)=(sum(Para.P(s_,:).*(squeeze((xePrimeMat(xind,Rind,:)))'-[x x]).^2)-EDeltaXMat(xind,Rind).^2)^0.5;
        EDeltaRMat(xind,Rind)=sum(Para.P(s_,:).*squeeze(RprimeMat(xind,Rind,:))')-R;
        VDeltaRMat(xind,Rind)=(sum(Para.P(s_,:).*(squeeze( (RprimeMat(xind,Rind,:)))'-[R R]).^2)-EDeltaRMat(xind,Rind).^2)^.5;
        xTarget=[x R s_];
    VDirect=funeval(c(s_,:)',V(s_),xTarget(1:2));
    CheckFit(xind,Rind)=(VDirect-V_new)/V_new;
    
ctr=ctr+1;        
    end
    



end
x0=SimulationData.InitialState_x;
R0=SimulationData.InitialState_R;
xSS=SimulationData.SteadyState_x;
RSS=SimulationData.SteadyState_R;



figure()
    [C,h] = contour(xPlotGrid,RPlotGrid,TransMat');
    xlabel('x')
    ylabel('$\rho$','Interpreter','Latex')
    text_handle = clabel(C,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on
    x1=[x0,xSS];
    y1=[R0,RSS];
    plot(x1(1),y1(1),'.','MarkerSize',30,'Color','k');
    hold on
    plot(x1(2),y1(2),'.','MarkerSize',30,'Color','r');
    hold on
colorbar('peer',gca);

fr=SimulationData.SampleFrequency;
T=length(SimulationData.xHist);
xHist=SimulationData.xHist;
RHist=SimulationData.RHist;
plot(xHist(1:fr:T),RHist(1:fr:T),'ok')

title('Transfers')
print('-dpng',[plotpath 'TransfersDynamics.png'])

figure()
[C,h] = contour(xPlotGrid,RPlotGrid,TauMat');
xlabel('x')
    ylabel('$\rho$','Interpreter','Latex')
    
    text_handle = clabel(C,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on
    x1=[x0,xSS];
    y1=[R0,RSS];
    plot(x1(1),y1(1),'.','MarkerSize',30,'Color','k');
    hold on
    plot(x1(2),y1(2),'.','MarkerSize',30,'Color','r');
    hold on
    
colorbar('peer',gca);
fr=SimulationData.SampleFrequency;
T=length(SimulationData.xHist);
xHist=SimulationData.xHist;
RHist=SimulationData.RHist;
plot(xHist(1:fr:T),RHist(1:fr:T),'ok')

title('Taxes')


print('-dpng',[plotpath 'TaxesDynamics.png'])

figure()
[C,h] = contour(xPlotGrid,RPlotGrid,BtildePrimeMat');
xlabel('x')
    ylabel('$\rho$','Interpreter','Latex')
    
    text_handle = clabel(C,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on
    x1=[x0,xSS];
    y1=[R0,RSS];
    plot(x1(1),y1(1),'.','MarkerSize',30,'Color','k');
    hold on
    plot(x1(2),y1(2),'.','MarkerSize',30,'Color','r');
    hold on
colorbar('peer',gca);

fr=SimulationData.SampleFrequency;
T=length(SimulationData.xHist);
xHist=SimulationData.xHist;
RHist=SimulationData.RHist;
plot(xHist(1:fr:T),RHist(1:fr:T),'ok')
title('Debt')

print('-dpng',[plotpath 'DebtDynamics.png'])

