function iter=MainBellman(Para)
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
[ Para,V] = BuidGrid( Para);
disp('Msg: Completed definition of functional space')

%% INITIALIZE THE COEFF
disp('Msg: Initializing the Value function....')

tic
[ x_state, c, PolicyRulesStore] = InitializeCoeff( Para, V)    ;
toc
disp('Msg: .... Completed')

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
u2btild_slice=x_state(:,1) ;
R_slice=x_state(:,2) ;
s_slice=x_state(:,3) ;
GridSize=Para.GridSize;
for iter=2:Para.Niter
    tic    
    % Clear the records for the arrays that store the index of the point in
    % the doman (w.r.t x_state) where the inner optimization failed
    IndxSolved=[];
    IndxUnSolved=[];
    ExitFlag=[];  
    
    % Initialize the initial guess for the policy rules that the inner
    % optimization will solve
    PolicyRulesStoreOld=PolicyRulesStore;
    parfor ctr=1:GridSize/2        
        u2btild=u2btild_slice(ctr) ;
        R=R_slice(ctr) ;
        s_=s_slice(ctr);
        % INITAL GUESS FOR THE INNER OPTIMIZATION
        xInit=PolicyRulesStore(ctr,:);
        % INNER OPTIMIZATION
        [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,xInit',Para,0);
        ExitFlag(ctr)=exitflag;
        VNew(ctr)=V_new;
        %UODATE POLICY RULES
        if exitflag==1
            PolicyRulesStore(ctr,:)=PolicyRules;
        end
        
    end
    % -- IID CASE -----
   % In the IID case we solve it for s=1 and use the solution to populate
   % s=2.
    ExitFlag(GridSize/2+1:GridSize)=ExitFlag(1:GridSize/2);
    VNew(GridSize/2+1:GridSize)=VNew(1:GridSize/2);
    PolicyRulesStore(GridSize/2+1:GridSize,:)=PolicyRulesStore(1:GridSize/2,:);
    % --- ----------------------
    
    % UNRESOLVED POINTS
    HandleUnresovledPoints
    % UPDATE NEW COEFF
    UpdateCoeffecients
    if length(IndxUnSolved)./GridSize >.01
        disp('exiting for a new grid')
        break;
    end
    
end
% STORE THE FINAL COEFF
save([ Para.datapath Para.StoreFileName] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','x_state','Para','V');
