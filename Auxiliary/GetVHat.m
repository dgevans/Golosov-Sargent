function [c,V] =GetVHat(Para,RGrid,c0SS,VSS,x_state,PolicyRulesStoreOld)
close all;



%% Build Grid for the state variables
% This setups up the functional space and the grid.
%u2btildMin=-(Para.theta_1-Para.theta_2)/(1-Para.beta)*(1/(Para.n1*Para.theta_1+Para.n2*Para.theta_2-Para.g(1)));
%u2btildMin=u2btildMin/4;
u2btildMin=Para.u2btildMin;
u2btildMax=Para.u2btildMax;
u2btildGrid=linspace(u2btildMin,u2btildMax,Para.u2btildGridSize);
Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;
RMin=RGrid.RMin;
RMax=RGrid.RMax;
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
   
    c=c0SS;
    
    % slicing the state space for parfor loop later
    u2btild_slice=x_state(:,1) ;
    R_slice=x_state(:,2) ;
    s_slice=x_state(:,3) ;
    % This stores the values of the policy functions and multipliers that last
    % worked
    
    %% ITERATE ON THE VALUE FUNCTION
    
    
    %Para.g
        tic
       
        IndxSolved=[];
        IndxUnSolved=[];
        ExitFlag=[];
        for ctr=1:GridSize/2
            
            u2btild=u2btild_slice(ctr) ;
            R=R_slice(ctr) ;
            s_=s_slice(ctr);
            [xInit]=GetInitialApproxPolicy([u2btild R s_],x_state,PolicyRulesStoreOld);
            %xInit=PolicyRulesStore(ctr,:);
            [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,VSS,xInit',Para,0);
            ExitFlag(ctr)=exitflag;
            VNew(ctr)=V_new;
            % update policyrules guess only for exitflag==1
            if exitflag==1
            PolicyRulesStore(ctr,:)=PolicyRules(1:3);
            else
                disp([u2btild_slice(ctr) R_slice(ctr) ])
                PolicyRulesStore(ctr,:)=xInit;
            end
        end
        
        ExitFlag(GridSize/2+1:GridSize)=ExitFlag(1:GridSize/2);
        VNew(GridSize/2+1:GridSize)=VNew(1:GridSize/2);
        PolicyRulesStore(GridSize/2+1:GridSize,:)=PolicyRulesStore(1:GridSize/2,:);         
        IndxUnSolved=find(~(ExitFlag==1));
        IndxSolved=find(ExitFlag==1);
        sprintf('fraction of nodes unresolved at the first pass = %1.3f',length(IndxUnSolved)./GridSize)
       
%         % --  Rsolve the FOC at points that failed in the first round -----
             NumTrials=5;
          sprintf('Resolving the unresolved points using alterative routine ')
 
             UnResolvedPoints
%             if NumResolved>0
%                 Numtrials=10;
%                 sprintf('Resolving the unresolved points using alterative routine ')
%                 UnResolvedPoints;
%             end
        IndxUnSolved=find(~(ExitFlag==1));
        IndxSolved=find(ExitFlag==1);
        IndxSolved_1=IndxSolved(IndxSolved<=GridSize/Para.sSize);
        IndxSolved_2=IndxSolved(IndxSolved>GridSize/Para.sSize);
        %
        % Obtain the new coeffecins by projecting the Cheb polynomials for
        % both the value functions
        
        
      cNew(1,:)=funfitxy(V(1),x_state(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
       
        
        cNew(2,:)=cNew(1,:);
     
        
        
        cOld=c;
        % update the guess by taking a weighted average of the old and new
        % coeffecients
        c=cNew*Para.grelax+(1-Para.grelax)*cOld;
        
        disp('Completed Iteration No - ')
  % Store the difference
        cdiff(1,:)=sum(abs(c-cNew))';
      
        toc
       x_state(IndxUnSolved,:)
save([ Para.datapath 'cVHat.mat'] , 'c','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','x_state','Para','V','cdiff');
