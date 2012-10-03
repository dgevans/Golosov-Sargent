function  PlotSimulationCommonshockAlt( X,T,SimTitle,K,gHist,plotpath,texpath)
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
Moments(i,K)=corr(X.data(startIndex:end,i),n1*theta_1Hist(startIndex:end,i)+n2*theta_2Hist(startIndex:end,i));
end


rowLabels = SimTitle;
columnLabels = {'Mean','Std','AutoCorr','Corr with $n_1theta_1+n_2theta_2$'};
matrix2latex(Moments, [texpath X.name 'Moments.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.Kf', 'size', 'tiny');



end
