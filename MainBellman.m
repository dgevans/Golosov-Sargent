function iter=MainBellman(Para,RGrid,InitData)
close all;
% This is the main file for computing the minimally stochastic case for BGP
% preferences
%% NOTATION
% x= u_2 btild
% R = u_2/u_1
% - -- - - - -
clc
close all
%%
%  This sets up the flag which controls how the bellman function is
%  initialized. if set to no - it uses the data passed by the used in
%  InitData

switch nargin
    case 0 % No arguments
        SetParaStruc % set Para
        flagComputeInitCoeff='yes';
        flagSetRGrid='no';
    case 1 % given Para
        flagComputeInitCoeff='yes';
        flagSetRGrid='no';
    case 2 % given para and RGrid
        flagComputeInitCoeff='yes';
        flagSetRGrid='yes';
    case 3 % given para , RGrid and data for initialization
        flagComputeInitCoeff='no';
        flagSetRGrid='yes';
        cInit=InitData.c;
        VInit=InitData.V;
end


%% Compute the undistorted FB
s_=1;


%% Build Grid for the state variables
% This setups up the functional space and the grid.
%u2btildMin=-(Para.theta_1-Para.theta_2)/(1-Para.beta)*(1/(Para.n1*Para.theta_1+Para.n2*Para.theta_2-Para.g(1)));
%u2btildMin=u2btildMin/4;
if isfield(Para,'flagSetu2BtildGrid')
    flagSetu2BtildGrid=Para.flagSetu2BtildGrid;
else
    flagSetu2BtildGrid=0;
end

u2btildMin=-2.5;
u2btildMax=2.5;
if flagSetu2BtildGrid==1
    disp('using user defined grid on x')
u2btildMin=Para.u2btildMin;
u2btildMax=Para.u2btildMax;
end

u2btildGrid=linspace(u2btildMin,u2btildMax,Para.u2btildGridSize);

Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;
Rbar0=1;
for u2btild_ind=1:Para.u2btildGridSize
    findR=@(R) getValueC1(u2btildGrid(u2btild_ind),R,s_,Para);
    Rbar(u2btild_ind)=fzero(findR,Rbar0);
    Rbar0=Rbar(u2btild_ind);
end


% R=u_2/u_1 = c1/c2


RMin=max(Rbar)*1.01;
RMax=max(Rbar)*1.7;
if strcmpi(flagSetRGrid,'yes')==1
    disp('setting RGrid with user inputs')
    RMin=RGrid.RMin;
    RMax=RGrid.RMax;
end
RGrid=linspace(RMin,RMax,Para.RGridSize);

Para.RGrid=RGrid;
GridSize=Para.u2btildGridSize*Para.RGridSize*Para.sSize;
Para.GridSize=GridSize;
Para.u2btildMin=u2btildMin;
Para.u2btildMax=u2btildMax;
Para.RMax=RMax;
Para.RMin=RMin;
%% Define the funtional space
V(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
V(2) = V(1);
GridPoints=[Para.u2btildLL Para.u2btildUL;RMin RMax];
rowLabels = {'$x$','$R$'};
columnLabels = {'Lower Bound','Upper Bounds'};
matrix2latex(GridPoints, [Para.texpath 'GridPoints.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.2f', 'size', 'tiny');

%% Set the Parallel Config
err=[];
try
    matlabpool('size')
catch err
end
if isempty(err)
    
    
    if(matlabpool('size') == 0)
        %matlabpool close
        matlabpool open local;
        
    end
    
    
end

%% INITIALIZE THE COEFF
%  This function computes c1,c2,l1,l2 and the value for an arbitrary x, R.
% This section solves for V i.e the value function at the end of period 1
% with g_t=g for all t >1. since the value function is static we need to
% solve a equation in c_1 for each x,R. Th function getValueC1 does the job
% by solving for the two roots of this equation and using the one that
% supports the highest utility
tic
gTrue=Para.g;
Para.g=mean(gTrue)*ones(2,1);
for s_=1:Para.sSize
    n=1;
    if s_==1
        
        
        for u2btildctr=1:Para.u2btildGridSize
            for Rctr=1:Para.RGridSize
                
                u2btild_=u2btildGrid(u2btildctr);
                R_=RGrid(Rctr);
                %if R_>Rbar(u2btildctr)
                x_state_(s_,n,:)=[u2btild_ R_ ];
                % Solve for  c1
                cRat = R_^(-1/Para.sigma);
                c1_1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(1))/(Para.n1+cRat*Para.n2);
                c1_2 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g(2))/(Para.n1+cRat*Para.n2);
                c2_1 = cRat*c1_1;
                options = optimset('Display','off');
                [xSS,~,exitFlag] = fsolve(@(x) SteadyStateResiduals(x,u2btild_,R_,Para,s_),[c1_1 c1_2 c2_1],options);
                [res, c1_, c2_, l1_, l2_] = SteadyStateResiduals(xSS,u2btild_,R_,Para,s_);
                if(exitFlag ~= 1)
                    R_
                    u2btild_
                    res
                end
                
                V0(s_,n) = (Para.alpha_1*uAlt(c1_,l1_,Para.psi,Para.sigma)+Para.alpha_2*uAlt(c2_,l2_,Para.psi,Para.sigma))*Para.P(s_,:)'/(1-Para.beta);
                
                
                xInit_0(s_,n,:)=[c1_(s_) c2_(s_) l1_(s_) l2_(s_) u2btild_ R_ u2btild_];
                n=n+1;
                %end
            end
        end
        c0(s_,:)=funfitxy(V(s_),squeeze(x_state_(s_,:,:)),squeeze(V0(s_,:))' );
    else
        c0(s_,:)=c0(s_-1,:);
        V0(s_,:)=V0(s_-1,:);
        xInit_0(s_,:)=xInit_0(s_-1,:);
        if strcmpi(flagComputeInitCoeff,'no')
            PolicyRulesStore(GridSize/2+1:GridSize,:)=PolicyRulesStore(1:GridSize/2,:);
        end
        
    end
end
%disp('Number of points solved in initialization')
%sum(ExitFlagT)
%disp('Number of points solved out of a total of ')
%length(ExitFlagT)

Para.g=gTrue;
x_state=vertcat([squeeze(x_state_(1,:,:)) ones(length(x_state_),1)] ,[squeeze(x_state_(1,:,:)) 2*ones(length(x_state_),1)]);
%scatter(x_state(:,1),x_state(:,2))
c=c0;
save([ Para.datapath 'c1.mat' ] , 'c');

% slicing the state space for parfor loop later
u2btild_slice=x_state(:,1) ;
R_slice=x_state(:,2) ;
s_slice=x_state(:,3) ;
% This stores the values of the policy functions and multipliers that last
% worked

PolicyRulesWorked=[xInit_0(1,1,1) xInit_0(2,1,1) xInit_0(1,1,2)];

% This stores the policy rules for each point in the state
% space.
PolicyRulesStore1=[squeeze(xInit_0(1,:,1))' squeeze(xInit_0(1,:,1))' ...
    squeeze(xInit_0(1,:,2))' squeeze(xInit_0(1,:,2))'...
    squeeze(xInit_0(1,:,3))' squeeze(xInit_0(1,:,3))' ...
    squeeze(xInit_0(1,:,4))' squeeze(xInit_0(1,:,4))' ....
    squeeze(xInit_0(1,:,5))' squeeze(xInit_0(1,:,5))' ....
    squeeze(xInit_0(1,:,6))' squeeze(xInit_0(1,:,6))' ....
    squeeze(xInit_0(1,:,7))' squeeze(xInit_0(1,:,7))' ....
    ];
PolicyRulesStore2=[squeeze(xInit_0(2,:,1))' squeeze(xInit_0(2,:,1))' ...
    squeeze(xInit_0(2,:,2))' squeeze(xInit_0(2,:,2))'...
    squeeze(xInit_0(2,:,3))' squeeze(xInit_0(2,:,3))' ...
    squeeze(xInit_0(2,:,4))' squeeze(xInit_0(2,:,4))' ....
    squeeze(xInit_0(2,:,5))' squeeze(xInit_0(2,:,5))' ....
    squeeze(xInit_0(2,:,6))' squeeze(xInit_0(2,:,6))' ....
    squeeze(xInit_0(2,:,7))' squeeze(xInit_0(2,:,7))' ....
    ];

if strcmpi(flagComputeInitCoeff,'yes')
    PolicyRulesStore=vertcat(PolicyRulesStore1,PolicyRulesStore2);
end
%% ITERATE ON THE VALUE FUNCTION


%Para.g
for iter=2:Para.Niter
    tic
    
    IndxSolved=[];
    IndxUnSolved=[];
    ExitFlag=[];
    PolicyRulesStoreOld=PolicyRulesStore;
    parfor ctr=1:GridSize/2
        
        u2btild=u2btild_slice(ctr) ;
        R=R_slice(ctr) ;
        s_=s_slice(ctr);
        %[xInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStoreOld);
        xInit=PolicyRulesStore(ctr,:);
        [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,xInit',Para,0);
        ExitFlag(ctr)=exitflag;
        VNew(ctr)=V_new;
        % update policyrules guess only for exitflag==1
        if exitflag==1
            PolicyRulesStore(ctr,:)=PolicyRules;
        end
        
    end
    
    ExitFlag(GridSize/2+1:GridSize)=ExitFlag(1:GridSize/2);
    VNew(GridSize/2+1:GridSize)=VNew(1:GridSize/2);
    PolicyRulesStore(GridSize/2+1:GridSize,:)=PolicyRulesStore(1:GridSize/2,:);
    IndxUnSolved=find(~(ExitFlag==1));
    IndxSolved=find(ExitFlag==1);
    sprintf('fraction of nodes unresolved at the first pass = %1.3f',length(IndxUnSolved)./GridSize)
    
    % --  Rsolve the FOC at points that failed in the first round -----
    if mod(iter,Para.ResolveCtr)==0
        NumTrials=5;
        disp('Points that failed the first round of FOC')
        x_state(IndxUnSolved,:)
        sprintf('Resolving the unresolved points using alterative routine ')
        UnResolvedPoints
        if NumResolved>0
            Numtrials=10;
            sprintf('Resolving the unresolved points using alterative routine ')
            UnResolvedPoints;
        end
    end
    IndxUnSolved=find(~(ExitFlag==1));
    IndxSolved=find(ExitFlag==1);
    IndxSolved_1=IndxSolved(IndxSolved<=GridSize/Para.sSize);
    IndxSolved_2=IndxSolved(IndxSolved>GridSize/Para.sSize);
    %
    % Obtain the new coeffecins by projecting the Cheb polynomials for
    % both the value functions
    
    
    cNew(1,:)=funfitxy(V(1),x_state(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
    
    
    cNew(2,:)=cNew(1,:);
    
    
    % Store the difference
    cdiff(iter,:)=sum(abs(c-cNew))';
    cOld=c;
    % update the guess by taking a weighted average of the old and new
    % coeffecients
    c=cNew*Para.grelax+(1-Para.grelax)*cOld;
    
    % Error in sup Norm
    ErrorInSupNorm(iter-1)=max(abs(VNew(IndxSolved_1)'-funeval(cOld(1,:)',V(1),x_state(IndxSolved_1,1:2))));
    
    disp('Completed Iteration No - ')
    disp(iter)
    toc
    
    if mod(iter,1)==0
        save([ Para.datapath  'c_' num2str(iter) '.mat' ] , 'c','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','x_state','Para','V');
    end
    if length(IndxUnSolved)./GridSize >.01
        disp('exiting for a new grid')
        break;
    end
    
end

save([ Para.datapath Para.StoreFileName] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','x_state','Para','V');
