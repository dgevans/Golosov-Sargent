
clear all


s_=1;

btild_1=-2

%% Calibrated Persistence Shocks
% GShocks
load ~/Golosov-Sargent/Data/Draft/cGShocks.mat
%PolicyFunctions
[ x0,R0] = solveTime0Problem(Para,c,V,btild_1,s_) ;NormDecompReal(1,:)=ComputeConditionalChangeBudgetConstraint(Para,c,V,s_,domain,PolicyRulesStore,x0,R0);

% TFP
load ~/Golosov-Sargent/Data/Draft/cTFP.mat
%PolicyFunctions
[ x0,R0] = solveTime0Problem(Para,c,V,btild_1,s_) ;NormDecompReal(2,:)=ComputeConditionalChangeBudgetConstraint(Para,c,V,s_,domain,PolicyRulesStore,x0,R0);

% TFPIneq
load ~/Golosov-Sargent/Data/Draft/cTFPIneq.mat
%PolicyFunctions
[ x0,R0] = solveTime0Problem(Para,c,V,btild_1,s_) ;NormDecompReal(3,:)=ComputeConditionalChangeBudgetConstraint(Para,c,V,s_,domain,PolicyRulesStore,x0,R0);

load ~/Golosov-Sargent/Data/Draft/cTFPIneqBetaShocks.mat
%PolicyFunctions
[ x0,R0] = solveTime0Problem(Para,c,V,btild_1,s_) ;NormDecompReal(4,:)=ComputeConditionalChangeBudgetConstraint(Para,c,V,s_,domain,PolicyRulesStore,x0,R0);


load ~/Golosov-Sargent/Data/Draft/cTFPIneqLargerBetaShocks.mat
%PolicyFunctions
[ x0,R0] = solveTime0Problem(Para,c,V,btild_1,s_) ;NormDecompReal(5,:)=ComputeConditionalChangeBudgetConstraint(Para,c,V,s_,domain,PolicyRulesStore,x0,R0);



 rowLabels = {'GShocks','TFP','TFP+Ineq','TFP+Ineq+Beta','TFP+Ineq+(Large) Beta'};
 columnLabels = {'$\Delta g$','$\Delta B$','$\Delta T$','$\Delta [\tau\theta_1l_1]$' ,'$\Delta [\tau\theta_2l_2]$' , '$\Delta Y$' ,'$\Delta \tau$'};
 matrix2latex( [NormDecompReal], 'DeompositionReal.tex' , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'normalsize');

