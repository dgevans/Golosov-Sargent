function PlotParallelSimulationsCommonShocks(SimDataPath,SimTexPath,SimPlotPath,SimTitle)

% This script plots the long simulation using data stored in
% Data/SimDataParallel.mat file. 


load( SimDataPath);

mkdir(SimPlotPath)
plotpath=SimPlotPath;
mkdir(SimTexPath)
texpath=SimTexPath;


K=size(u2btildHist,2);
T=100;

%SimTitle = {'Benchmark Simulation','Different Pareto Weights Simulation','Mean Preserving Spread in g Simulation', 'Increased Inequality Simulation'};




% -- u2btild ----------------------------------------------------------
X.data=u2btildHist;
X.sHist=sHist;
X.ylabel='x';
X.name ='x';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)




% -- R ----------------------------------------------------------
X.data=RHist;
X.sHist=sHist;
X.ylabel='R';
X.name ='R';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)


% -- FrishElaticity ----------------------------------------------------------
X.data=(Para.n1*(1./l1Hist-1)+Para.n2*(1./l2Hist-1))/(Para.n1+Para.n2);
X.sHist=sHist;
X.ylabel='FE';
X.name ='FrishElaticity';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)

% -- GDP ----------------------------------------------------------
X.data=YHist;
X.sHist=sHist;
X.ylabel='Y';
X.name ='Y';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)


% -- labor taxes ----------------------------------------------------------
X.data=TauHist;
X.sHist=sHist;
X.ylabel='tau';
X.name ='LaborTaxes';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)
   

% -- btild ----------------------------------------------------------
X.data=btildHist;
X.sHist=sHist;
X.ylabel='b2~';
X.name ='RelativeAssetsAgent2';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)


% -- Trans ----------------------------------------------------------
X.data=TransHist;
X.sHist=sHist;
X.ylabel='T';
X.name ='Transfers';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)


% -- AfterTaxIncomeAgent1 ----------------------------------------------------------
X.data=AfterTaxWageIncome_Agent1Hist;
X.sHist=sHist;
X.ylabel='After-tax wage income';
X.name ='AfterTaxWageIncomeAgent1';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)

% -- AfterTaxIncomeAgent2 ----------------------------------------------------------
X.data=AfterTaxWageIncome_Agent2Hist;
X.sHist=sHist;
X.ylabel='After-tax wage income';
X.name ='AfterTaxWageIncomeAgent2';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist,theta_2Hist,plotpath,Para.n1,Para.n2,texpath)

% -- IncomeFromAssetsAgent1 ----------------------------------------------------------
X.data=IncomeFromAssets_Agent1Hist;
X.sHist=sHist;
X.ylabel='Asset Income';
X.name ='IncomeFromAssetsAgent1';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)


% -- RtBt ----------------------------------------------------------
X.data=-IncomeFromAssets_Agent1Hist+btildHist(1:end-1,:);
X.sHist=sHist;
X.ylabel='RtBt';
X.name ='RtBt';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)

% -- Int Rates ----------------------------------------------------------
X.data=IntHist;
X.sHist=sHist;
X.ylabel='Int';
X.name ='Int';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)


% -- Gini Coeff ----------------------------------------------------------
X.data=GiniCoeffHist(2:end,:);
X.sHist=sHist(2:end,:);
X.ylabel='Gini';
X.name ='Gini';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(2:end,:),theta_2Hist(2:end,:),plotpath,Para.n1,Para.n2,texpath)

% -- Trans Diff----------------------------------------------------------
X.data=TransDiffHist;
X.sHist=sHist;
X.ylabel='TransDiff';
X.name ='TransDiff';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)


% -- LaborTaxAgent1Diff Diff----------------------------------------------------------
X.data=LaborTaxAgent1DiffHist;
X.sHist=sHist;
X.ylabel='LaborTaxAgent1Diff';
X.name ='LaborTaxAgent1Diff';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)

% -- LaborTaxAgent2Diff Diff----------------------------------------------------------
X.data=LaborTaxAgent2DiffHist;
X.sHist=sHist;
X.ylabel='LaborTaxAgent2Diff';
X.name ='LaborTaxAgent2Diff';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)

% -- DebtDiff Diff----------------------------------------------------------
X.data=DebtDiffHist;
X.sHist=sHist;
X.ylabel='DebtDiff';
X.name ='DebtDiffDiff';
PlotSimulationCommonshockAlt( X,T,SimTitle,K,theta_1Hist(1:end-1,:),theta_2Hist(1:end-1,:),plotpath,Para.n1,Para.n2,texpath)





%Government financing 
%--- Government Financing ------------------------------------------------
a=5000;
b=10000;

for i = 1:K
   
shock=sHist(a+1:b,i);
sDiff=sHist(a+1:b,i)-sHist(a:b-1,i);
ChangeNewAssets=(Para.n1*btildHist(a+1:b,i))-Para.n1*btildHist(a:b-1,i);
ChangeOldAssets=(Para.n1*btildHist(a:b-1,i))-Para.n1*btildHist(a-1:b-2,i);
ChangeAssetsDiff=(Para.n1*btildHist(a+1:b,i)-Para.n1*btildHist(a:b-1,i))-(Para.n1*btildHist(a:b-1,i)-Para.n1*btildHist(a-1:b-2,i));
TransDiff=(Para.n1+Para.n2)*(TransHist(a+1:b,i)-TransHist(a:b-1,i));
LaborTaxRevenueAgent1Diff=TauHist(a+1:b,i).*(Para.n1*theta_1Hist(a+1:b,i).*l1Hist(a+1:b,i))- TauHist(a:b-1,i).*(Para.n1*theta_1Hist(a+1:b,i).*l1Hist(a:b-1,i));
LaborTaxRevenueAgent2Diff=TauHist(a+1:b,i).*(Para.n2*theta_2Hist(a+1:b,i).*l2Hist(a+1:b,i))- TauHist(a:b-1,i).*(Para.n2*theta_2Hist(a+1:b,i).*l2Hist(a:b-1,i));
NetIntOnAssets=(Para.n1*btildHist(a:b-1,i).*(IntHist(a:b-1,i)-1))-(Para.n1*btildHist(a-1:b-2,i).*(IntHist(a-1:b-2,i)-1));
GovFundingComponents=[sDiff ChangeNewAssets ChangeOldAssets ChangeAssetsDiff TransDiff LaborTaxRevenueAgent1Diff LaborTaxRevenueAgent2Diff NetIntOnAssets];
IndxIncrease=find(GovFundingComponents(:,1)>0);
IndxDecrease=find(GovFundingComponents(:,1)<0);
IndxNoChange=find(GovFundingComponents(:,1)==0);

IndxNoChangeHigh=find(and(shock==2, GovFundingComponents(:,1)==0));
IndxNoChangeLow=find(and(shock==1, GovFundingComponents(:,1)==0));
GovFinTableChange=[mean(GovFundingComponents(IndxIncrease,1:end));mean(GovFundingComponents(IndxDecrease,1:end))]*100;
GovFinTableNoChange=[mean(GovFundingComponents(IndxNoChangeLow,1:end));mean(GovFundingComponents(IndxNoChangeHigh,1:end))]*100;
CheckGovFin=GovFinTableChange(:,1)+GovFinTableChange(:,2)+GovFinTableChange(:,3)-GovFinTableChange(:,4)-GovFinTableChange(:,5)-GovFinTableChange(:,6);
Decomp=[DebtDiffHist(a:b-1,i) 2*TransDiffHist(a:b-1,i) LaborTaxAgent1DiffHist(a:b-1,i) LaborTaxAgent2DiffHist(a:b-1,i)];
IndxGLow=find(sHist(a:b-1,i)==1) ;
DecompLow=mean(Decomp(IndxGLow,:)); 
IndxGHigh=find(sHist(a:b,i)==2) ;
DecompHigh=mean(Decomp(IndxGHigh,:)); 
 
rowLabels = {'$s\_=1$','$s\_=2$'};
columnLabels = {'$n_1[b''(h)-b''(l)]$','$(n_1+n_2)[T(h)-T(l)]$','$n_1\theta_1[l_1(h)\tau(h)-l_1(l)\tau(l)]$', '$n_2\theta_2[l_2(h)\tau(h)-l_2(l)\tau(l)]$'};
matrix2latex([ DecompLow;DecompHigh], [texpath 'GovHighLow' num2str(i) '.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');



rowLabels = {'(s=1->s=2)',' (s=2->s=1)'};
columnLabels = {'$\Delta s$','$\Delta b_{t+1}$','$\Delta b_{t}$','$\Delta (b_{t+1}-b_{t})$','$\Delta T$','$\Delta (\tau n_1\theta_1 l_1 )$','$\Delta (\tau n_2\theta_2 l_2)$','$\Delta ([\mathcal{R}-1]b_t)$'};
matrix2latex(GovFinTableChange, [texpath 'GovFinChange' num2str(i) '.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');

rowLabels = {'No Change (low)','No Change (high)'};
columnLabels = {'$\Delta s$','$\Delta b_{t+1}$','$\Delta b_{t}$','$\Delta (b_{t+1}-b_t)$','$\Delta T$','$\Delta (\tau n_1\theta_1 l_1 )$','$\Delta (\tau n_2\theta_2 l_2)$','$\Delta ([\mathcal{R}-1]b_t)$'};
matrix2latex(GovFinTableNoChange, [texpath 'GovFinNoChange' num2str(i) '.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');




end


end