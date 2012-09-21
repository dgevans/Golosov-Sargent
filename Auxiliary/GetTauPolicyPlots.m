% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function GetTauPolicyPlots(u2bdiffFineGrid,R,s_,iter,Para,plotpath)
close all;
olddatapath=Para.datapath;
oldtexpath=Para.texpath;
datapath=olddatapath;
load([datapath 'c' num2str(iter) '.mat'])
disp('Govt Exp')
g=Para.g

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

 
 for u2btildctr=1:length(u2bdiffFineGrid)
 
     u2btild=u2bdiffFineGrid(u2btildctr);
   [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag,fvec]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
    if exitflag==1
        IndxPrint(u2btildctr)=1;
    else
        IndxPrint(u2btildctr)=0;
    end
    c1=PolicyRules(1:2);
 c2=PolicyRules(3:4);
 l1=PolicyRules(5:6);
 l2=PolicyRules(7:8);
 ul2=(1-psi)./(1-l2);
uc2=psi./c2;

% TAU - From the WAGE optimality of Agent 2
Tau(u2btildctr,:)=1-(ul2./(theta_2.*uc2));
Trans(u2btildctr,:)=c2-l2.*ul2./uc2;
u2BtildePrime(u2btildctr,:)=PolicyRules(end-1:end);
 end
 
 
 
 
 figure()

 plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),1)- u2bdiffFineGrid(logical(IndxPrint))','k','LineWidth',2)
 hold on
 plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),2)-u2bdiffFineGrid(logical(IndxPrint))',':k','LineWidth',2)
    legend('g_l','g_h')
 xlabel('$x$','Interpreter','Latex')
 ylabel('$x(s)-x$','Interpreter','Latex')
 title(['$R=$' num2str(R)])

 
print(gcf,'-dpng',[plotpath 'Delta_X.png'])
 
 
 figure()
 subplot(1,2,1)
  plot(u2bdiffFineGrid(logical(IndxPrint)), Tau(logical(IndxPrint),1),'k','LineWidth',2)
 hold on
 plot(u2bdiffFineGrid(logical(IndxPrint)), Tau(logical(IndxPrint),2),':k','LineWidth',2)
 

     legend('g_l','g_h')

 xlabel('$x$','Interpreter','Latex')
 ylabel('$\tau$','Interpreter','Latex')
 title(['$R=$' num2str(R)])

 
subplot(1,2,2)
  plot(u2bdiffFineGrid(logical(IndxPrint)), Trans(logical(IndxPrint),1),'k','LineWidth',2)
 hold on
 plot(u2bdiffFineGrid(logical(IndxPrint)), Trans(logical(IndxPrint),2),':k','LineWidth',2)
 

     legend('g_l','g_h')

 xlabel('$x$','Interpreter','Latex')
 ylabel('$\tau$','Interpreter','Latex')
 title(['$R=$' num2str(R)]) 
 print(gcf,'-dpng',[plotpath 'Taxes_Transfers_Policy.png'])
end
 