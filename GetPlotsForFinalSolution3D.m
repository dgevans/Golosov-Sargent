% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetPlotsForFinalSolution3D(Para,Domain)
close all
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

surf(xPlotGrid,RPlotGrid,Values','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
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

% SOLVE THE T-0 PROBLEM given btild(-1)
c10guess=1;
c20guess=.5;
btild_1=Para.btild_1;
s_=1;

disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,s_,Para,c,V),[ c10guess c20guess],options);
if ~(exitflagv0==1)
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,s_,Para,c,V),[ 1 1/Para.RMax],options);
end

if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
    opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
        'GradObj','off','GradConstr','off',...
        'MaxIter',1000, ...
        'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
    lb=[0.001 0.001];
    ub=[10 10];
    [x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,s_,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
    %[x,~,exitflagv0,output,lambda]  =ktrlink(@(x) getValue0(x, btild_1,1,Para,c,V),[ c10guess c20guess],[],[],[],[],lb,ub,[],opts);
    
end
c10 = x(1);
c20 = x(2);
R0=(c10/c20)^(sigma);
TotalResources=(c10*n1+c20*n2+g(s_));
DenL2=theta_2*R0*n1+theta_2*n2;
l20=(TotalResources-theta_1*n1+ theta_2*n1*R0)/(DenL2);
l10= 1-(1-l20)*theta_2/theta_1*R0;
u2btildprime0=-(c20-c10)*(psi*c20^(-sigma))-((l10/(1-l10))*R0-l20/(1-l20))*(1-psi)+btild_1*psi*c20^(-sigma);
Rprime0=c20^(-sigma)/c10^(-sigma);



options = optimset('Display','off','TolX',1e-10);

 xState=fsolve(@(x) GetCrossingPoints(x,1,c,V,PolicyRulesStore,x_state,Para),[0 4]);
u2btild_=xState(1);
R_=xState(2);



figure()
subplot(3,1,1)
    [C,h] = contour(xPlotGrid,RPlotGrid,squeeze(TransMat(:,:,1))');
    text_handle = clabel(C,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on
    x1=[u2btildprime0,u2btild_]
    y1=[Rprime0,R_]
    plot(x1(1),y1(1),'.','MarkerSize',20,'Color','k');
    hold on
    plot(x1(2),y1(2),'.','MarkerSize',20,'Color','r');
    hold on
colorbar('peer',gca);

plot(x1,y1,':r','LineWidth',1)
title('Transfers')

subplot(3,1,2)
    [C,h] = contour(xPlotGrid,RPlotGrid,squeeze(TauMat(:,:,1))');
    text_handle = clabel(C,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on
    x1=[u2btildprime0,u2btild_]
    y1=[Rprime0,R_]
    plot(x1(1),y1(1),'.','MarkerSize',20,'Color','k');
    hold on
    plot(x1(2),y1(2),'.','MarkerSize',20,'Color','r');
    hold on
colorbar('peer',gca);

plot(x1,y1,':r','LineWidth',1)

title('Taxes')


subplot(3,1,3)
    [C,h] = contour(xPlotGrid,RPlotGrid,squeeze(BtildePrimeMat(:,:,1))');
    text_handle = clabel(C,h);
    set(text_handle,'BackgroundColor',[1 1 .6],...
        'Edgecolor',[.7 .7 .7])
    hold on
    x1=[u2btildprime0,u2btild_]
    y1=[Rprime0,R_]
    plot(x1(1),y1(1),'.','MarkerSize',20,'Color','k');
    hold on
    plot(x1(2),y1(2),'.','MarkerSize',20,'Color','r');
    hold on
colorbar('peer',gca);

plot(x1,y1,':r','LineWidth',1)

title('Debt')


figure()
surf(xPlotGrid,RPlotGrid,squeeze(TauMat(:,:,1))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
hold on
surf(xPlotGrid,RPlotGrid,squeeze(TauMat(:,:,2))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
axis tight
xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'Tau3D.png'])



figure()
surf(xPlotGrid,RPlotGrid,squeeze(TransMat(:,:,1))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
hold on
surf(xPlotGrid,RPlotGrid,squeeze(TransMat(:,:,2))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

xlabel('x')
ylabel('R')
print(gcf,'-dpng',[plotpath 'Tau3D.png'])


figure()
surf(xPlotGrid,RPlotGrid,squeeze(Deltau2BtildePrimeMat(:,:,1))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
hold on
surf(xPlotGrid,RPlotGrid,squeeze(Deltau2BtildePrimeMat(:,:,2))','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')

xlabel('x')
ylabel('R')
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
ylabel('R')
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
ylabel('R')
print(gcf,'-dpng',[plotpath 'VDeltaX3D.png'])

figure()
surf(xPlotGrid,RPlotGrid,VDeltaRMat','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
xlabel('x')
ylabel('R')
ylabel('R')
print(gcf,'-dpng',[plotpath 'VDeltaR3D.png'])


figure()
surf(xPlotGrid,RPlotGrid,CheckFit','FaceColor','interp',...
 'EdgeColor','none',...
 'FaceLighting','phong')
xlabel('x')
ylabel('R')
zlabel('FitError')
print(gcf,'-dpng',[plotpath 'ChebError3D.png'])


figure()
plot(ErrorInSupNorm)
xlabel('Iteration');
ylabel('Sup Norm Difference');
print(gcf,'-dpng',[plotpath 'CoeffConvergenceSupNorm.png'])

