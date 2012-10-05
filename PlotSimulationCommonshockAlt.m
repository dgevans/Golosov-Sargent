function PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,n1,n2,texpath)



BurnSampleRatio=.5;% Percentage of simulations to disregard
%%
%Long Simulations
figure()
for i = 1:K
    subplot(K,1,i)
    plot(X.data(:,i))
    xlabel('t')
    ylabel(X.ylabel,'Interpreter','Latex')
    title([X.name ' - Long Run Plot of ' SimTitle{i}]);
end
print(gcf,'-depsc2 ',[plotpath 'LongSimulations' X.name '.eps'])
print(gcf,'-dpng ',[plotpath 'LongSimulations' X.name '.png'])

%%
%Short Simulation

figure()
for i = 1:K
    subplot(K,1,i)
    XX.Data=X.data(2:T+1,i);
    XX.sHist=X.sHist(2:T+1,i);
    XX.name=X.ylabel;  
    PlotSimul(XX,1);
    title([X.name ' - First 100 periods of' SimTitle{i}],'Interpreter','Latex');
end

print(gcf,'-depsc2 ',[plotpath 'TruncSimulations' X.name 'First100.eps'])
print(gcf,'-dpng ',[plotpath 'TruncSimulations' X.name 'First100.png'])


%Short Simulation
figure()
for i = 1:K
    subplot(K,1,i)
    XX.Data=X.data(end-T+1:end,i);
    XX.sHist=X.sHist(end-T+1:end,i);
    XX.name=X.ylabel;  
    PlotSimul(XX,1);
    title([X.name ' - Last 100 periods of' SimTitle{i}])
end

print(gcf,'-depsc2 ',[plotpath 'TruncSimulations' X.name 'Last100.eps'])
print(gcf,'-dpng ',[plotpath 'TruncSimulations' X.name 'Last100.png'])


%%
% -- moments -------------------------------------------------------------
for i = 1:K
    startIndex = floor(BurnSampleRatio*length(X.data(:,i)));
Moments(i,1) =mean(X.data(startIndex:end,i));
Moments(i,2)=std(X.data(startIndex:end,i));
Moments(i,3)=corr(X.data(startIndex:end,i),X.data(startIndex-1:end-1,i));
Moments(i,4)=corr(X.data(startIndex:end,i),theta_1Hist(startIndex:end,i));
end


rowLabels = SimTitle;
columnLabels = {'Mean','Std','AutoCorr','Corr'};
matrix2latex(Moments, [texpath X.name 'Moments.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');



end