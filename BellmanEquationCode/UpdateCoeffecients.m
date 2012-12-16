% UPDATE THE NEW COEFFs
   
    cNew(1,:)=funfitxy(V(1),domain(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
    cNew(2,:)=cNew(1,:);
        
   % Store the difference
    cdiff(iter,:)=sum(abs(c-cNew))';
    cOld=c;
    % update the guess by taking a weighted average of the old and new
    % coeffecients
    c=cNew*Para.grelax+(1-Para.grelax)*cOld;
    
    % Error in sup Norm
    ErrorInSupNorm(iter-1)=max(abs(VNew(IndxSolved_1)'-funeval(cOld(1,:)',V(1),domain(IndxSolved_1,1:2))));
    disp('Completed Iteration No - ')
    disp(iter)
    toc
    
    if mod(iter,1)==0
        save([ Para.datapath  'c_' num2str(iter) '.mat' ] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','domain','Para','V');
    end
    