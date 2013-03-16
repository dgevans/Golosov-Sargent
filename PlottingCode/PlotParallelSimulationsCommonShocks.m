function PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle,Para)

% This script plots the long simulation using data stored in
% Data/SimDataParallel.mat file. 

BigT=Para.BigT;
% BigT controls the length of sample for long simulations plots

SimData=load( SimDataPath);
SimData=SimData.SimData;
K=length(SimData);
for k=1:K
        sHist(:,k)=[SimData(k).sHist];
        if ~(length(Para.g)>1)
Theta_1Hist(:,k)= SimData(k).Theta_1Hist;
Theta_2Hist(:,k)= SimData(k).Theta_2Hist;
ThetaShockDiffHist(:,k)= SimData(k).ThetaShockDiffHist;
        else
            gHist(:,k)=SimData(k).gHist;
            GShockDiffHist(:,k)= SimData(k).GShockDiffHist;
        end
        
xHist(:,k)= SimData(k).xHist;
RHist(:,k)= SimData(k).RHist;
TauHist(:,k)= SimData(k).TauHist;
YHist(:,k)= SimData(k).YHist;
TransHist(:,k)= SimData(k).TransHist;
btildHist(:,k)= SimData(k).btildHist;
c1Hist(:,k)= SimData(k).c1Hist;
c2Hist(:,k)= SimData(k).c2Hist;
l1Hist(:,k)= SimData(k).l1Hist;
l2Hist(:,k)= SimData(k).l2Hist;
IntHist(:,k)= SimData(k).IntHist;
IncomeFromAssets_Agent1Hist(:,k)= SimData(k).IncomeFromAssets_Agent1Hist;
AfterTaxWageIncome_Agent1Hist(:,k)= SimData(k).AfterTaxWageIncome_Agent1Hist;
AfterTaxWageIncome_Agent2Hist(:,k)= SimData(k).AfterTaxWageIncome_Agent2Hist;

TransDiffHist(:,k)= SimData(k).TransDiffHist;
LaborTaxAgent1DiffHist(:,k)= SimData(k).LaborTaxAgent1DiffHist;
LaborTaxAgent2DiffHist(:,k)= SimData(k).LaborTaxAgent2DiffHist;
DebtDiffHist(:,k)= SimData(k).DebtDiffHist;
GiniCoeffHist(:,k)= SimData(k).GiniCoeffHist;


end
mkdir(SimPlotPath)
plotpath=SimPlotPath;
mkdir(SimTexPath)
texpath=SimTexPath;



T=100;

%SimTitle = {'Benchmark Simulation','Different Pareto Weights Simulation','Mean Preserving Spread in g Simulation', 'Increased Inequality Simulation'};




% -- x ----------------------------------------------------------
X.data=xHist(2:end,:);
X.sHist=sHist;
X.ylabel='x';
X.name ='x';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)



% -- R ----------------------------------------------------------
X.data=RHist;
X.sHist=sHist;
X.ylabel='R';
X.name ='R';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)






% -- GDP ----------------------------------------------------------
X.data=YHist;
X.sHist=sHist;
X.ylabel='Y';
X.name ='Y';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)

% -- g/y ----------------------------------------------------------
if ~(length(Para.g)>1)
X.data=Para.g./YHist;
else
    X.data=gHist./YHist;
end
X.sHist=sHist;
X.ylabel='g/y';
X.name ='gyratio';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)



% -- labor taxes ----------------------------------------------------------
X.data=TauHist;
X.sHist=sHist;
X.ylabel='tau';
X.name ='LaborTaxes';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)
   

% -- btild ----------------------------------------------------------
X.data=btildHist;
X.sHist=sHist;
X.ylabel='b2~';
X.name ='RelativeAssetsAgent2';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)


% -- Trans ----------------------------------------------------------
X.data=TransHist;
X.sHist=sHist;
X.ylabel='T';
X.name ='Transfers';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)


% -- AfterTaxIncomeAgent1 ----------------------------------------------------------
X.data=AfterTaxWageIncome_Agent1Hist;
X.sHist=sHist;
X.ylabel='After-tax wage income';
X.name ='AfterTaxWageIncomeAgent1';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)

% -- AfterTaxIncomeAgent2 ----------------------------------------------------------
X.data=AfterTaxWageIncome_Agent2Hist;
X.sHist=sHist;
X.ylabel='After-tax wage income';
X.name ='AfterTaxWageIncomeAgent2';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)

% -- IncomeFromAssetsAgent1 ----------------------------------------------------------
X.data=IncomeFromAssets_Agent1Hist;
X.sHist=sHist;
X.ylabel='Asset Income';
X.name ='IncomeFromAssetsAgent1';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)


% -- RtBt ----------------------------------------------------------
X.data=-IncomeFromAssets_Agent1Hist+btildHist(1:end-1,:);
X.sHist=sHist;
X.ylabel='RtBt';
X.name ='RtBt';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)




% -- Int Rates ----------------------------------------------------------
X.data=IntHist;
X.sHist=sHist;
X.ylabel='Int';
X.name ='Int';
PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)


end
% -- FrishElaticity ----------------------------------------------------------
%X.data=(Para.n1*(1./l1Hist-1)+Para.n2*(1./l2Hist-1))/(Para.n1+Para.n2);
%X.sHist=sHist;
%X.ylabel='FE';
%X.name ='FrishElaticity';
%PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)

% % -- Gini Coeff ----------------------------------------------------------
% X.data=GiniCoeffHist;
% X.sHist=sHist;
% X.ylabel='Gini';
% X.name ='Gini';
% PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist,plotpath,texpath)
% 
% % -- Trans Diff----------------------------------------------------------
% X.data=TransDiffHist;
% X.sHist=sHist;
% X.ylabel='TransDiff';
% X.name ='TransDiff';
% PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)
% 
% 
% % -- LaborTaxAgent1Diff Diff----------------------------------------------------------
% X.data=LaborTaxAgent1DiffHist;
% X.sHist=sHist;
% X.ylabel='LaborTaxAgent1Diff';
% X.name ='LaborTaxAgent1Diff';
% PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)
% 
% % -- LaborTaxAgent2Diff Diff----------------------------------------------------------
% X.data=LaborTaxAgent2DiffHist;
% X.sHist=sHist;
% X.ylabel='LaborTaxAgent2Diff';
% X.name ='LaborTaxAgent2Diff';
% PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)
% 
% % -- DebtDiff Diff----------------------------------------------------------
% X.data=DebtDiffHist;
% X.sHist=sHist;
% X.ylabel='DebtDiff';
% X.name ='DebtDiffDiff';
% PlotSimulationCommonshockAlt( X,T,BigT,SimTitle,K,YHist(1:end-1,:),plotpath,texpath)
% 
% 

%
% % Government financing 
% % --- Government Financing ------------------------------------------------
% a=length(RHist)-5000;
% b=length(RHist);
% %a=1000
% %b=2000
% for i = 1:K
%    
%     gshock=sHist(a+1:b,i);
% gDiff=sHist(a+1:b,i)-sHist(a:b-1,i);
% ChangeNewAssets=(Para.n1*btildHist(a+1:b,i))-Para.n1*btildHist(a:b-1,i);
% ChangeOldAssets=(Para.n1*btildHist(a:b-1,i))-Para.n1*btildHist(a-1:b-2,i);
% ChangeAssetsDiff=(Para.n1*btildHist(a+1:b,i)-Para.n1*btildHist(a:b-1,i))-(Para.n1*btildHist(a:b-1,i)-Para.n1*btildHist(a-1:b-2,i));
% TransDiff=(Para.n1+Para.n2)*(TransHist(a+1:b,i)-TransHist(a:b-1,i));
% LaborTaxRevenueAgent1Diff=TauHist(a+1:b,i).*(Para.n1*Para.theta_1*l1Hist(a+1:b,i))- TauHist(a:b-1,i).*(Para.n1*Para.theta_1*l1Hist(a:b-1,i));
% LaborTaxRevenueAgent2Diff=TauHist(a+1:b,i).*(Para.n2*Para.theta_2*l2Hist(a+1:b,i))- TauHist(a:b-1,i).*(Para.n2*Para.theta_2*l2Hist(a:b-1,i));
% NetIntOnAssets=(Para.n1*btildHist(a:b-1,i).*(IntHist(a:b-1,i)-1))-(Para.n1*btildHist(a-1:b-2,i).*(IntHist(a-1:b-2,i)-1));
% GovFundingComponents=[gDiff ChangeNewAssets ChangeOldAssets ChangeAssetsDiff TransDiff LaborTaxRevenueAgent1Diff LaborTaxRevenueAgent2Diff NetIntOnAssets];
% IndxIncrease=find(GovFundingComponents(:,1)>0);
% IndxDecrease=find(GovFundingComponents(:,1)<0);
% IndxNoChange=find(GovFundingComponents(:,1)==0);
% 
% ghigh=max(sHist(:,i));
% glow=min(sHist(:,i));
% 
% IndxNoChangeHigh=find(and(gshock==ghigh, GovFundingComponents(:,1)==0));
% IndxNoChangeLow=find(and(gshock==glow, GovFundingComponents(:,1)==0));
% GovFinTableChange=[mean(GovFundingComponents(IndxIncrease,1:end));mean(GovFundingComponents(IndxDecrease,1:end))]*100;
% GovFinTableNoChange=[mean(GovFundingComponents(IndxNoChangeLow,1:end));mean(GovFundingComponents(IndxNoChangeHigh,1:end))]*100;
% CheckGovFin=GovFinTableChange(:,1)+GovFinTableChange(:,2)+GovFinTableChange(:,3)-GovFinTableChange(:,4)-GovFinTableChange(:,5)-GovFinTableChange(:,6);
% Decomp=[ GShockDiffHist(a:b-1,i) DebtDiffHist(a:b-1,i) 2*TransDiffHist(a:b-1,i) LaborTaxAgent1DiffHist(a:b-1,i) LaborTaxAgent2DiffHist(a:b-1,i)];
% IndxGLow=find(sHist(a:b-1,i)==glow) ;
% DecompLow=mean(Decomp(IndxGLow,:)); 
% IndxGHigh=find(sHist(a:b-1,i)==ghigh) ;
% DecompHigh=mean(Decomp(IndxGHigh,:)); 
%  
% rowLabels = {'$g\_=g_l$','$g\_=g_h$'};
% columnLabels = {'$g_h-g_l$','$n_1[b''(h)-b''(l)]$','$(n_1+n_2)[T(h)-T(l)]$','$n_1\theta_1[l_1(h)\tau(h)-l_1(l)\tau(l)]$', '$n_2\theta_2[l_2(h)\tau(h)-l_2(l)\tau(l)]$'};
% matrix2latex([ DecompLow;DecompHigh], [texpath 'GovHighLow' num2str(i) '.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');
% 
% 
% 
% rowLabels = {'Increase (low->high)','Decrease (high->low)'};
% columnLabels = {'$\Delta g$','$\Delta b_{t+1}$','$\Delta b_{t}$','$\Delta (b_{t+1}-b_{t})$','$\Delta T$','$\Delta (\tau n_1\theta_1 l_1 )$','$\Delta (\tau n_2\theta_2 l_2)$','$\Delta ([\mathcal{R}-1]b_t)$'};
% matrix2latex(GovFinTableChange, [texpath 'GovFinChange' num2str(i) '.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');
% 
% rowLabels = {'No Change (low)','No Change (high)'};
% columnLabels = {'$\Delta g$','$\Delta b_{t+1}$','$\Delta b_{t}$','$\Delta (b_{t+1}-b_t)$','$\Delta T$','$\Delta (\tau n_1\theta_1 l_1 )$','$\Delta (\tau n_2\theta_2 l_2)$','$\Delta ([\mathcal{R}-1]b_t)$'};
% matrix2latex(GovFinTableNoChange, [texpath 'GovFinNoChange' num2str(i) '.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');
% 
% 
% 
% 
% end
