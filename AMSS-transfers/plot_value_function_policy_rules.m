 
    load('~/Golosov-Sargent/Data/temp/AMSS_Solution_tr_2shocks.mat')

    psi=Para.psi;
    beta=Para.beta;
    g=Para.g;
    pi=Para.pi
    xGrid=Para.xGrid
    plotpath='/Graphs'
    
    % Figure 1 : Value function
    figure()
    s_=1;
    plot(xGrid,funeval(coeff(s_,:)',V(s_),xGrid),':k','LineWidth',2)
    hold on
    s_=2;
    plot(xGrid,funeval(coeff(s_,:)',V(s_),xGrid),'k','LineWidth',2)
    xlabel('x')
    ylabel('V(x)')
    
    print(gcf,'-dpng',[ plotpath 'FigAMSSValueFunction1.png'])

    % Figure 2a : Policy Rule -n
    figure()
    s_=1;
    for s=1:Para.sSize
    plot(xGrid,squeeze(n(:,s_,s)),'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('n(x)')
    
    title('$s_{-1}=1$','Interpreter','Latex')

    print(gcf,'-dpng',[ plotpath 'FigAMSSLaborPolicy1.png'])

    % Figure 2a : Policy Rule - taxes
    % tau=1-ul/uc = 1- (1-psi)/()
    
    tax=@(n,s) 1-((1-psi)/(psi)).*(n-g(s))./(1-n);
    figure()
    s_=1;
    for s=1:sSize
    plot(xGrid,tax(squeeze(n(:,s_,s)),s),'LineWidth',2)
    hold on
    end
    xlabel('x')
    ylabel('\tau(x)')
        
    title('$s_{-1}=1$','Interpreter','Latex')

    print(gcf,'-dpng',[ plotpath 'FigAMSSTaxPolicy.png'])


    % Figure 2b : Policy Rule -x'

     figure()
     s_=1;
     for s=1:sSize
     plot(xGrid,squeeze(xprime(:,s_,s))-xGrid,'LineWidth',2)
     
     hold on
     end
    
     print(gcf,'-dpng',[ plotpath 'FigAMSSxPolicy1.png'])
    

    