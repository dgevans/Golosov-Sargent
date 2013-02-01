%{ In order to call this you need to be in the following position:

% execute RunMainWithAltSigmas and break on line 75 in MainBellman.
% set ctr = 1 and copy/paste 77-81 of MainBellman into command window
% set breakpoint on line 40 in CheckGradNag
% copy/paste CheckGradNAG(x,R,s_,c,V,xInit',Para) into command window
% Set break point somewhere in BelObjectiveUncondGradNAGBGP.
% call this function (copy/paste temp_slyon into command window)

%}

BelObjectiveUncondGradNAGBGP(0, zInit, 0, 0)