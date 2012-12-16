IndxUnSolved=find(~(ExitFlag==1));
    IndxSolved=find(ExitFlag==1);
  %disp('Total Unresolved Points')
  sprintf(' Unresolved so far  %1.2f',length(find(IndxUnSolved)))
NumUnsolved=length(IndxUnSolved);
for i=1:NumUnsolved
    IndxSolved=find(ExitFlag==1);
        uns_indx=IndxUnSolved(i);
    
     x=x_slice(uns_indx) ;
        R=R_slice(uns_indx) ;
        s_=s_slice(uns_indx);
    disp('Resolving...');
       
disp([x R s_]);

        
        
%% TRY I
NumTrials;
x0=[];
[PolicyRulesInit,xref]=GetInitialApproxPolicy([x R,s_],domain(IndxSolved,:),PolicyRulesStore(IndxSolved,:));
 x0(1,:)=linspace(xref(1),x,NumTrials) ;
 x0(2,:)=linspace(xref(2),R,NumTrials) ;
for tr_indx=1:NumTrials
 [PolicyRules, V_new,exitflag,~]=CheckGradNAG(x0(1,tr_indx),x0(2,tr_indx),s_,c,V,PolicyRulesInit,Para);        
 PolicyRulesInit=PolicyRules;
end




if ~(exitflag==1)
%% TRY 2
xguesss2=[x R,s_].*(1+(sign([x R,s_]-xref)*.05));
x0=[];
[PolicyRulesInit,xref]=GetInitialApproxPolicy(xguesss2,domain(IndxSolved,:),PolicyRulesStore(IndxSolved,:));
 x0(1,:)=linspace(xref(1),x,NumTrials) ;
 x0(2,:)=linspace(xref(2),R,NumTrials) ;
for tr_indx=1:NumTrials
 [PolicyRules, V_new,exitflag,~]=CheckGradNAG(x0(1,tr_indx),x0(2,tr_indx),s_,c,V,PolicyRulesInit,Para);        
 PolicyRulesInit=PolicyRules;
end
end

 
 
        ExitFlag(uns_indx)=exitflag;
        VNew(uns_indx)=V_new;
        if exitflag==1
        PolicyRulesStore(uns_indx,:)=PolicyRules;
        end
           IndxSolved=find(ExitFlag==1);

        
end
           

    IndxSolved=find(ExitFlag==1);
    IndxUnSolved=find(~(ExitFlag==1));
   disp('Number of points resolved with alternative guess')
   NumResolved=NumUnsolved-length(IndxUnSolved)
   
 