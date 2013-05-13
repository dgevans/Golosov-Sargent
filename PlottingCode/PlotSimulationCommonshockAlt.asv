function  PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)
BurnSampleRatio=.5;% Percentage of simulations to disregard

%%
%Long Simulations
figure()
for i = 1:K
    subplot(K,1,i)
    plot(X.data(2:BigT,i),'k','LineWidth',2)
    xlabel('t')
    ylabel(X.ylabel,'Interpreter','Latex')
    title([X.name ' - Long Run Plot  ' SimTitle{i}]);
end
print(gcf,'-depsc2 ',[plotpath 'LongSimulations' X.name '.eps'])
print(gcf,'-dpng ',[plotpath 'LongSimulations' X.name '.png'])

%%
%Short Simulation
%Short Simulation

figure()
for i = 1:K
    subplot(K,1,i)
    XX.Data=X.data(2:T+1,i);
    XX.sHist=X.sHist(2:T+1,i);
    XX.name=X.ylabel;  
    PlotSimul(XX,1);
    title([X.name ' - First 100 periods ' SimTitle{i}],'Interpreter','Latex');
end

print(gcf,'-depsc2 ',[plotpath 'TruncSimulations' X.name 'First100.eps'])
print(gcf,'-dpng ',[plotpath 'TruncSimulations' X.name 'First100.png'])


figure()
for i = 1:K
    subplot(K,1,i)
    XX.Data=X.data(end-T+1:end,i);
    XX.sHist=X.sHist(end-T+1:end,i);
    XX.name=X.ylabel;  
    PlotSimul(XX,1);
    title([X.name ' - Last 100 periods ' SimTitle{i}])
end

print(gcf,'-depsc2 ',[plotpath 'TruncSimulations' X.name 'Last100.eps'])
print(gcf,'-dpng ',[plotpath 'TruncSimulations' X.name 'Last100.png'])

% 
% %%
 % -- moments -------------------------------------------------------------
 % last TT periods
 TT=100;
 for i = 1:K
   
 MomentsLastTT(i,1) =mean(X.data(end-TT:end,i));
 MomentsLastTT(i,2)=std(X.data(end-TT:end,i));
 MomentsLastTT(i,3)=corr(X.data(end-TT:end,i),X.data(end-TT-1:end-1,i));
 MomentsLastTT(i,4)=corr(X.data(end-TT:end,i),YHist(end-TT:end,i));
 end
 % First TT periods
 t1=10;
 for i = 1:K
     
 MomentsFirstTT(i,1) =mean(X.data(t1:t1+TT,i));
 MomentsFirstTT(i,2)=std(X.data(t1:t1+TT,i));
 MomentsFirstTT(i,3)=corr(X.data(t1:t1+TT,i),X.data(t1-1:t1+TT-1,i));
 MomentsFirstTT(i,4)=corr(X.data(t1:t1+TT,i),YHist(t1:t1+TT,i));
 end

 
% 
% 
 rowLabels = {'First 100 periods','Last 100 periods'};
 columnLabels = {'Mean','Std','AutoCorr','Corr with g'};
 matrix2latex( [MomentsFirstTT;MomentsLastTT], [texpath X.name 'Moments.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');
 

end
