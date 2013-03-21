% This script documents the disgnostic checks for the FOC

close all
clear all


% ---------------------------------------------------------------------------------------
%% CHECK THE FOCs
% ---------------------------------------------------------------------------------------


% LOAD THE COEFF
LastIter=250;
load(['Data/c' num2str(LastIter) '.mat'])
SetParaStruc
GetPlots(2,LastIter,Para)
% Pick up the test points
NumTestPoints=20;
s_=1;
% FIX the controls at FB levels
[c1FB c2FB l1FB l2FB yFB g_yFB_h Agent1WageShareFB_h]=getFB(Para,1);
TestControls=[c1FB c1FB c2FB];

% delta h for the numerical differences
deltah=.00001; 

u2btildLL=Para.u2btildMin*.9;
u2btildUL=Para.u2btildMax*.9;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];

% Compute the DIFF in the numerical and analytic derivatives

for n=1:NumTestPoints
    u2btild=ucbtild_bounds(1)+(ucbtild_bounds(2)-ucbtild_bounds(1))*rand;
R=Rbounds(1)+(Rbounds(2)-Rbounds(1))*rand;     
xTarget(n,:)=[u2btild R s_];
DiffDerivaties(n,:)=CheckDerivatives(u2btild,R,s_,c,V,TestControls,Para,deltah);
end
matrix2latex([xTarget], [Para.texpath 'TargetPoints.tex'] , 'rowLabels', {}, 'columnLabels', {'x','R','$s\_$'}, 'alignment', 'c', 'format', '%-6.2g', 'size', 'tiny');
matrix2latex([DiffDerivaties], [Para.texpath 'CheckDerivatives.tex'] , 'rowLabels', {}, 'columnLabels', {'$c_1(1)$','$c_1(2)$','$c_2(1)$'}, 'alignment', 'c', 'format', '%-6.2g', 'size', 'tiny');






% Summarize the finding in a TEX table



% ---------------------------------------------------------------------------------------
%% Check the Const optimization routine
% ---------------------------------------------------------------------------------------
% Compute the zero of ComputeFOCUnc and ComputeFOCCons
for n=1:NumTestPoints
    u2btild=xTarget(n,1);
R=xTarget(n,2);
s_=xTarget(n,3);
[PolicyRulesInit]=GetInitialApproxPolicy(xTarget(n,:),x_state,PolicyRulesStore);
 [DiffFOCZero(n,:),DiffFOCOpt(n,:),exitflag(n,:)]=CheckFOC(u2btild,R,s_,c,V,PolicyRulesInit,Para) ;
end
matrix2latex([DiffFOCZero], [Para.texpath 'CheckFOCCons.tex'] , 'rowLabels', {}, 'columnLabels', {'$c_1(1)$','$c_1(2)$','$c_2(1)$'}, 'alignment', 'c', 'format', '%-6.2g', 'size', 'tiny');
matrix2latex([DiffFOCOpt], [Para.texpath 'CheckFOCwithOptmizer.tex'] , 'rowLabels', {}, 'columnLabels', {'$c_1(1)$','$c_1(2)$','$c_2(1)$'}, 'alignment', 'c', 'format', '%-6.2g', 'size', 'tiny');
matrix2latex([exitflag], [Para.texpath 'ExitFlag.tex'] , 'rowLabels', {}, 'columnLabels', {'ExitFlag'}, 'alignment', 'c', 'format', '%-6.2g', 'size', 'tiny');
% Summarize the difference in a TEx table


% ---------------------------------------------------------------------------------------
%Solve and Store the results for the deterministic case
% ---------------------------------------------------------------------------------------

% Change the paramters to the determisitic case
SetParaStruc
gTrue=Para.g;
Para.g=mean(gTrue)*ones(2,1);
Para.datapath=[Para.datapath 'Deterministic/'];
Para.texpath=[Para.texpath 'Deterministic/'];
Para.plotpath=[Para.plotpath 'Deterministic/'];
mkdir(Para.datapath);
mkdir(Para.plotpath);
mkdir(Para.texpath);
Para.Niter=5;
StartIter=1;
MainBellman(Para,StartIter)
GetPlots(2,Para.Niter,Para)


% ---------------------------------------------------------------------------------------
%FOC Errors
% ---------------------------------------------------------------------------------------

% Load the coeff from the last iteration
% Plot the FOC residuals




% ---------------------------------------------------------------------------------------
%Compute the Policy Rules for the first iteration
% ---------------------------------------------------------------------------------------

% Load the coeff from the first iteration
LastIter=2;
load(['Data/c' num2str(LastIter) '.mat'])
SetParaStruc
Para.plotpath=[Para.plotpath 'FirstIter/'];
mkdir(Para.datapath);
mkdir(Para.plotpath);
mkdir(Para.texpath);
GetPlots(2,2,Para)
% Plot the Policy Rules 







