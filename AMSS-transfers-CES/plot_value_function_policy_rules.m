 close all
 clear all
 
load('~/Golosov-Sargent/Data/temp/AMSS_ces_tr_2_shocks.mat')
sol_tr_2=AMSS_ces_tr_2_shocks;
load('~/Golosov-Sargent/Data/temp/AMSS_ces_tr_3_shocks.mat')
sol_tr_3=AMSS_ces_tr_3_shocks;


  %  plotpath='/home/anmol/Dropbox/2011RA/FiscalPolicy/GolosovProjectCode/SteadyStateNotes/AMSS/Images/'
    C={':k','k','-.k'};
    % Figure 1 : Value function
    s_=1;
    figure()
    subplot(1,2,1)
    for s=1:sol_tr_2.Para.sSize
    plot(funeval(sol_tr_2.coeff(s_,:)',sol_tr_2.V(s_),sol_tr_2.xGrid(1:end-3),1),funeval(sol_tr_2.coeff(s,:)',sol_tr_2.V(s),sol_tr_2.xprime((1:end-3),s_,s),1),C{s},'LineWidth',2)
    hold on
    end
    title(' Transfers - 2 shocks')
    xlabel('x')
    ylabel('\mu(x)')
    
    subplot(1,2,2)
    for s=1:sol_tr_3.Para.sSize
    plot(funeval(sol_tr_3.coeff(s_,:)',sol_tr_3.V(s_),sol_tr_3.xGrid(1:end-3),1),funeval(sol_tr_3.coeff(s,:)',sol_tr_3.V(s),sol_tr_3.xprime((1:end-3),s_,s),1),C{s},'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('\mu(x)')
    
    title(' Transfers - 3 shocks')
    
    print(gcf,'-dpng',[ plotpath 'FigAMSSMultiplier_tr.png'])

    % Figure 2a : Policy Rule -n
    figure()
    subplot(1,2,1)
    s_=1;
    for s=1:sol_tr_2.Para.sSize
    plot(sol_tr_2.xGrid,squeeze(sol_tr_2.n(:,s_,s)),C{s},'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('n(x)')
  title(' Labor : Transfers - 2 shocks')
  legend('g_l','g_h')  
    subplot(1,2,2)
     s_=1;
    for s=1:sol_tr_3.Para.sSize
    plot(sol_tr_3.xGrid,squeeze(sol_tr_3.n(:,s_,s)),C{s},'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('n(x)')
   
  title(' Labor :Transfers - 3 shocks')
  legend('g_l','g_m','g_h')  
  
    print(gcf,'-dpng',[ plotpath 'FigAMSSLaborPolicy_tr.png'])

    

    % Figure 2b : Policy Rule -x'

     figure()
     subplot(1,2,1)
     s_=1;
     for s=1:sol_tr_2.Para.sSize
     plot(sol_tr_2.xGrid,squeeze(sol_tr_2.xprime(:,s_,s))-sol_tr_2.xGrid,C{s},'LineWidth',2)
     
     hold on
     end
     xlabel('x')
     ylabel('$\Delta x(s)$','Interpreter','Latex')
     title(' Transfers - 2 shocks')
  legend('g_l','g_h')
     subplot(1,2,2)
    s_=1;
     for s=1:sol_tr_3.Para.sSize
     plot(sol_tr_3.xGrid,squeeze(sol_tr_3.xprime(:,s_,s))-sol_tr_3.xGrid,C{s},'LineWidth',2)
     
     hold on
     end
     title(' Transfers - 3 shocks')
     xlabel('x')
     ylabel('$\Delta x(s)$','Interpreter','Latex')
     
     print(gcf,'-dpng',[ plotpath 'FigAMSSxPolicy_tr.png'])
     
     legend('g_l','g_m','g_h')
    figure()
    T=10000;
    subplot(1,2,1)
plot(sol_tr_2.SimData.xHist(2:T),'k','LineWidth',2)
title(' Transfers - 2 shocks')

subplot(1,2,2)
plot(sol_tr_3.SimData.xHist(2:T),'k','LineWidth',2)

title(' Transfers - 3 shocks')
print(gcf,'-dpng',[ plotpath 'FigAMSSSamplePath_tr.png'])



close all
 clear all
 
load('~/Golosov-Sargent/Data/temp/AMSS_ces_no_tr_2_shocks.mat')
sol_no_tr_2=AMSS_ces_no_tr_2_shocks;
load('~/Golosov-Sargent/Data/temp/AMSS_ces_no_tr_3_shocks.mat')
sol_no_tr_3=AMSS_ces_no_tr_3_shocks;


    plotpath='/home/anmol/Dropbox/2011RA/FiscalPolicy/GolosovProjectCode/SteadyStateNotes/AMSS/Images/'
    C={':k','k','-.k'};
    % Figure 1 : Value function
    s_=1;
    figure()
    subplot(1,2,1)
    for s=1:sol_no_tr_2.Para.sSize
    plot(funeval(sol_no_tr_2.coeff(s_,:)',sol_no_tr_2.V(s_),sol_no_tr_2.xGrid(1:end-3),1),funeval(sol_no_tr_2.coeff(s,:)',sol_no_tr_2.V(s),sol_no_tr_2.xprime((1:end-3),s_,s),1),C{s},'LineWidth',2)
    hold on
    end
    title(' No Transfers - 2 shocks')
    xlabel('x')
    ylabel('\mu(x)')
    subplot(1,2,2)
    for s=1:sol_no_tr_3.Para.sSize
    plot(funeval(sol_no_tr_3.coeff(s_,:)',sol_no_tr_3.V(s_),sol_no_tr_3.xGrid(1:end-3),1),funeval(sol_no_tr_3.coeff(s,:)',sol_no_tr_3.V(s),sol_no_tr_3.xprime((1:end-3),s_,s),1),C{s},'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('\mu(x)')
    
    title(' No Transfers - 3 shocks')
    
    print(gcf,'-dpng',[ plotpath 'FigAMSSMultiplier_no_tr.png'])
    
    
    % Figure 2a : Policy Rule -n
    figure()
    subplot(1,2,1)
    s_=1;
    for s=1:sol_no_tr_2.Para.sSize
    plot(sol_no_tr_2.xGrid,squeeze(sol_no_tr_2.n(:,s_,s)),C{s},'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('n(x)')
  title(' Labor : No Transfers - 2 shocks')
  legend('g_l','g_h')  
    subplot(1,2,2)
     s_=1;
    for s=1:sol_no_tr_3.Para.sSize
    plot(sol_no_tr_3.xGrid,squeeze(sol_no_tr_3.n(:,s_,s)),C{s},'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('n(x)')
   
  title(' Labor : No Transfers - 3 shocks')
  legend('g_l','g_m','g_h')  
  
    print(gcf,'-dpng',[ plotpath 'FigAMSSLaborPolicy_no_tr.png'])

    

    % Figure 2b : Policy Rule -x'

     figure()
     subplot(1,2,1)
     s_=1;
     for s=1:sol_no_tr_2.Para.sSize
     plot(sol_no_tr_2.xGrid,squeeze(sol_no_tr_2.xprime(:,s_,s))-sol_no_tr_2.xGrid,C{s},'LineWidth',2)
     
     hold on
     end
     xlabel('x')
     ylabel('$\Delta x(s)$','Interpreter','Latex')
     title(' No Transfers - 2 shocks')
  legend('g_l','g_h')
     subplot(1,2,2)
    s_=1;
     for s=1:sol_no_tr_3.Para.sSize
     plot(sol_no_tr_3.xGrid,squeeze(sol_no_tr_3.xprime(:,s_,s))-sol_no_tr_3.xGrid,C{s},'LineWidth',2)
     
     hold on
     end
     title(' No Transfers - 3 shocks')
     xlabel('x')
     ylabel('$\Delta x(s)$','Interpreter','Latex')
     
     print(gcf,'-dpng',[ plotpath 'FigAMSSxPolicy_no_tr.png'])
     
     legend('g_l','g_m','g_h')
    figure()
    T=10000;
    subplot(1,2,1)
plot(sol_no_tr_2.SimData.xHist(2:T),'k','LineWidth',2)
title(' No Transfers - 2 shocks')

subplot(1,2,2)
plot(sol_no_tr_3.SimData.xHist(2:T),'k','LineWidth',2)

title(' No Transfers - 3 shocks')

print(gcf,'-dpng',[ plotpath 'FigAMSSSamplePath_no_tr.png'])
