clear all
clc
close all
AMSSTransfers=load('~/Golosov-Sargent/Data/temp/AMSS_Solution_tr.mat');
AMSSNoTransfers=load('~/Golosov-Sargent/Data/temp/AMSS_Solution.mat');
AMSSNoTransfers.xGrid=AMSSNoTransfers.Para.xGrid;
    psi=AMSSTransfers.Para.psi;
    beta=AMSSTransfers.Para.beta;
    g=AMSSTransfers.Para.g;
    pi=AMSSTransfers.Para.pi;

% We now compute the grid for x. 
for s_=1:2
n_fb=g+psi*(1-g);
c_fb=psi*(1-g);
uc_fb=1./(1-g);
Euc_fb=sum(pi(s_,:).*uc_fb);
int_fb=uc_fb./(beta*Euc_fb);
x_fb(s_,:)=(((1-psi)/(1)).*((1).*n_fb./(1-n_fb))-psi).*(1-int_fb).^(-1);
end
Para.n_fb=n_fb;
Para.c_fb=c_fb;
Para.x_fb=x_fb(1,1);
    
   ssPol=[.7 .7 0 -0.9];
    get_root_ss_nag= @(num,ssPol,user,iflag) getSteadyState(num,ssPol,AMSSTransfers.Para,user,iflag)  ;
    [ssPol,~,exitflag]=c05qb(get_root_ss_nag,ssPol,'xtol',1e-10);
    n_ss=ssPol(1:2);
    phi_ss=ssPol(3);
    x_ss=ssPol(4);
    mkdir('Graphs')
    plotpath='Graphs/';
    
    % Figure 1 : Value function
    figure()
    s_=1;
    plot(AMSSTransfers.xGrid,funeval(AMSSTransfers.coeff(s_,:)',AMSSTransfers.V(s_),AMSSTransfers.xGrid),'k','LineWidth',2)
    hold on
    s_=2;
    plot(AMSSNoTransfers.xGrid,funeval(AMSSNoTransfers.coeff(s_,:)',AMSSNoTransfers.V(s_),AMSSNoTransfers.xGrid),'b','LineWidth',2)
    xlabel('x')
    ylabel('V(x)')
    vline([Para.x_fb x_ss],{'r',':r'})

    print(gcf,'-dpng',[ plotpath 'FigCompAMSSValueFunction.png'])

    % Figure 2 : Policy Rule -n
    figure()
    s_=1;
    subplot(1,2,1)
    plot(AMSSTransfers.xGrid,squeeze(AMSSTransfers.n(:,s_,1)),':k','LineWidth',2)
    hold on
    plot(AMSSTransfers.xGrid,squeeze(AMSSTransfers.n(:,s_,2)),'k','LineWidth',2)
    xlabel('x')
    ylabel('n(x)')
    vline([Para.x_fb x_ss],{'r',':r'})
    
    title('Transfers','Interpreter','Latex')

    subplot(1,2,2)
    s_=2;
    plot(AMSSNoTransfers.xGrid,squeeze(AMSSNoTransfers.n(:,s_,1)),':b','LineWidth',2)
    hold on
    plot(AMSSNoTransfers.xGrid,squeeze(AMSSNoTransfers.n(:,s_,2)),'b','LineWidth',2)
    xlabel('x')
    ylabel('n(x)')
    title('No Transfers','Interpreter','Latex')
        vline([Para.x_fb x_ss],{'r',':r'})

    print(gcf,'-dpng',[ plotpath 'FigCompAMSSLaborPolicy.png'])

        % Figure 2 : Policy Rule -c
    figure()
    s_=1;
    subplot(1,2,1)
    plot(AMSSTransfers.xGrid,squeeze(AMSSTransfers.n(:,s_,1))-AMSSTransfers.Para.g(1),':k','LineWidth',2)
    hold on
    plot(AMSSTransfers.xGrid,squeeze(AMSSTransfers.n(:,s_,2))-AMSSTransfers.Para.g(2),'k','LineWidth',2)
    xlabel('x')
    ylabel('c(x)')
    vline([Para.x_fb x_ss],{'r',':r'})
    
    title('Transfers','Interpreter','Latex')

    subplot(1,2,2)
    s_=2;
    plot(AMSSNoTransfers.xGrid,squeeze(AMSSNoTransfers.n(:,s_,1))-AMSSNoTransfers.Para.g(1),':b','LineWidth',2)
    hold on
    plot(AMSSNoTransfers.xGrid,squeeze(AMSSNoTransfers.n(:,s_,2))-AMSSNoTransfers.Para.g(2),'b','LineWidth',2)
    xlabel('x')
    ylabel('c(x)')
    title('No Transfers','Interpreter','Latex')
        vline([Para.x_fb x_ss],{'r',':r'})

    print(gcf,'-dpng',[ plotpath 'FigCompAMSSConsumptionPolicy.png'])

    
    % Figure 3 : Policy Rule - taxes
    % tau=1-ul/uc = 1- (1-psi)/()
    tax=@(n,s) 1-((1-psi)/(psi)).*(n-g(s))./(1-n);
    figure()
    s_=1;
    subplot(1,2,1)
    plot(AMSSTransfers.xGrid,tax(squeeze(AMSSTransfers.n(:,s_,1)),1),':k','LineWidth',2)
    hold on
    plot(AMSSTransfers.xGrid,tax(squeeze(AMSSTransfers.n(:,s_,2)),2),'k','LineWidth',2)
    xlabel('x')
    ylabel('\tau(x)')
        vline([Para.x_fb x_ss],{'r',':r'})

    title('Transfers','Interpreter','Latex')

    subplot(1,2,2)
    s_=2;
    plot(AMSSNoTransfers.xGrid,tax(squeeze(AMSSNoTransfers.n(:,s_,1)),1),':b','LineWidth',2)
    hold on
    plot(AMSSNoTransfers.xGrid,tax(squeeze(AMSSNoTransfers.n(:,s_,2)),2),'b','LineWidth',2)
    xlabel('x')
    ylabel('\tau(x)')
    title('No Transfers','Interpreter','Latex')
        vline([Para.x_fb x_ss],{'r',':r'})

    print(gcf,'-dpng',[ plotpath 'FigCompAMSSTaxPolicy.png'])


    % Figure 2b : Policy Rule -x'

     figure()
     s_=1;
     subplot(1,2,1)
     plot(AMSSTransfers.xGrid,squeeze(AMSSTransfers.xprime(:,s_,1))-AMSSTransfers.xGrid,':k','LineWidth',2)
     hold on
     plot(AMSSTransfers.xGrid,squeeze(AMSSTransfers.xprime(:,s_,2))-AMSSTransfers.xGrid,'k','LineWidth',2)
    vline([Para.x_fb x_ss],{'r',':r'})
title('Transfers')
 xlabel('x')
  ylabel('$\Delta x(s)$','Interpreter','Latex')
           
hold on

     s_=2;
     subplot(1,2,2)
     plot(AMSSNoTransfers.xGrid,squeeze(AMSSNoTransfers.xprime(:,s_,1))-AMSSNoTransfers.xGrid,':b','LineWidth',2)
     hold on
     plot(AMSSNoTransfers.xGrid,squeeze(AMSSNoTransfers.xprime(:,s_,2))-AMSSNoTransfers.xGrid,'b','LineWidth',2)
     hold on
     print(gcf,'-dpng',[ plotpath 'FigCompAMSSxPolicy.png'])
    vline([Para.x_fb x_ss],{'r',':r'})
xlabel('x')
 ylabel('$\Delta x(s)$','Interpreter','Latex')
 
title('No Transfers')

     
     
    % Figure 3 : Policy Rule :x'-T
%x'-x-Tuc = n./(1-n)*((1-psi)/psi)+x*(R-1)-psi

computeDeltaX_Tr =@(x,n) (n./(1-n)) *(1-psi)+x*(1./(beta*(n-g).*sum(pi(s_,:).*(1./(n-g)))))-x-psi;
for xind=1:length(AMSSTransfers.Para.xGrid)
    DeltaXTr(xind,:)=computeDeltaX_Tr(AMSSTransfers.xGrid(xind),squeeze(AMSSTransfers.n(xind,1,:))');

end
figure()
subplot(1,2,1)
plot(AMSSTransfers.xGrid,DeltaXTr(:,1),':k','LineWidth',2)
hold on
plot(AMSSTransfers.xGrid,DeltaXTr(:,2),'k','LineWidth',2)
vline(x_fb(1),':r')
    vline([Para.x_fb x_ss],{'r',':r'})
xlabel('x')
  ylabel('$\Delta x(s)-ucT$','Interpreter','Latex')
 
title('Transfers')


for xind=1:length(AMSSNoTransfers.xGrid)
    DeltaXTr(xind,:)=computeDeltaX_Tr(AMSSNoTransfers.xGrid(xind),squeeze(AMSSNoTransfers.n(xind,1,:))');

end
subplot(1,2,2)
plot(AMSSNoTransfers.xGrid,DeltaXTr(:,1),':b','LineWidth',2)
hold on
plot(AMSSNoTransfers.xGrid,DeltaXTr(:,2),'b','LineWidth',2)
vline(x_fb(1),':r')
 vline([Para.x_fb x_ss],{'r',':r'})
 title('No Transfers')
 xlabel('x')
 ylabel('$\Delta x(s)-ucT$','Interpreter','Latex')
 
   print(gcf,'-dpng',[ plotpath 'FigCompAMSSx_trPolicy.png'])


load('~/Golosov-Sargent/Data/temp/AMSSTransfers.mat')
load('~/Golosov-Sargent/Data/temp/AMSSNoTransfers.mat')
figure()
plot(AMSSNoTransfers.SimData.bzero.xHist,'b','LineWidth',2)
hold on
plot(AMSSTransfers.SimData.bzero.xHist,'k','LineWidth',2)
   print(gcf,'-dpng',[ plotpath 'FigCompAMSSSimulationbzero.png'])

figure()
plot(AMSSNoTransfers.SimData.balt.xHist,'b','LineWidth',2)
hold on
plot(AMSSTransfers.SimData.balt.xHist,'k','LineWidth',2)
print(gcf,'-dpng',[ plotpath 'FigCompAMSSSimulationbalt.png'])


