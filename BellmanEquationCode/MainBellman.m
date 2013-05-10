function iter=MainBellman(Para,BellmanData)
close all;
% This is the main file computes the value function via time iteration for
% the parameters passsed in the structure Para. 


%% NOTATION
% x= u_2 btild
% R = u_2/u_1
% - -- - - - -
clc
close all

% DEFAULT PARAMETERS
switch nargin
    case 0 % No arguments
        %% Set the Parallel Config
err=[];
try
    matlabpool('size')
catch err
end
disp('Msg: Using default parameters stored in SetParaStruc')
        SetParaStruc % set Para
end    

% BUILD GRID
[ Para,V,xhat,Rhat] = BuidGrid( Para);
disp('Msg: Completed definition of functional space')
%% INITIALIZE THE COEFF

if nargin==2
    
disp('Msg: Initializing the Value function with existing coeff')

tic
[ domain, c, PolicyRulesStore] = InitializeCoeffWithExistingCoeff( Para, V,BellmanData);
toc
else

disp('Msg: Initializing the Value function....')

tic
[ domain, c, PolicyRulesStore] = InitializeCoeff( Para, V)    ;
toc
disp('Msg: .... Completed')
end
%% OPEN MATLAB PARALLEL WORKERS
err=[];
try
    matlabpool('size')
catch err
end

if isempty(err)
   
    if(matlabpool('size') == 0)
        matlabpool open local;
        
    end
    
    
end

%% ITERATE ON THE VALUE FUNCTION
% This block iterates on the bellman equation

% Slicing the state space for parfor loop
x_slice=domain(:,1) ;
R_slice=domain(:,2) ;
s_slice=domain(:,3) ;
GridSize=Para.GridSize;
ErrorInSupNorm(1)=1;
S = length(Para.P);
for iter=2:Para.Niter
    tic    
    % Clear the records for the arrays that store the index of the point in
    % the doman (w.r.t domain) where the inner optimization failed
    IndxSolved=[];
    IndxUnSolved=[];
    ExitFlag=[];  
    
    % Initialize the initial guess for the policy rules that the inner
    % optimization will solve
    PolicyRulesStoreOld=PolicyRulesStore;
    %parfor ctr=1:GridSize       
    for ctr=1:GridSize/S
    %xInit=PolicyRulesStore(1,:);
    %for ctr=1:GridSize/S       
        x=x_slice(ctr) ;
        R=R_slice(ctr) ;
        s_=s_slice(ctr);
        % INITAL GUESS FOR THE INNER OPTIMIZATION
        xInit=PolicyRulesStore(ctr,:);
        % INNER OPTIMIZATION
        %[PolicyRules, V_new,exitflag,~]=CheckGradNAG2Shocks(x,R,s_,c,V,xInit',Para);
        [PolicyRules, V_new,exitflag,~]=CheckGradNAG(x,R,s_,c,V,xInit',Para);
        ExitFlag(ctr)=exitflag;
        VNew(ctr)=V_new;
        %UODATE POLICY RULES
        if exitflag==1
            PolicyRulesStore(ctr,:)=PolicyRules;
        end
        %xInit=PolicyRules;
    end
    % -- IID CASE -----
   % In the IID case we solve it for s=1 and use the solution to populate
   % s=2.
   for s=1:S-1
    ExitFlag(s*GridSize/S+1:(s+1)*GridSize/S)=ExitFlag(1:GridSize/S);
    VNew(s*GridSize/S+1:(s+1)*GridSize/S)=VNew(1:GridSize/S);
    PolicyRulesStore(s*GridSize/S+1:(s+1)*GridSize/S,:)=PolicyRulesStore(1:GridSize/S,:);
   end
    % --- ----------------------
    
    % UNRESOLVED POINTS
    HandleUnresovledPoints
    % UPDATE NEW COEFF
    UpdateCoeffecients
    if length(IndxUnSolved)./GridSize >.01
        disp('exiting for a new grid')
        break;
    end
    
    if ErrorInSupNorm(iter-1) < Para.ctol;
        disp('convergence criterion met')
        break;
    end
    
end

save([ Para.datapath Para.StoreFileName] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','domain','Para','V','xhat','Coeff_xhat','Rhat','Coeff_Rhat');
