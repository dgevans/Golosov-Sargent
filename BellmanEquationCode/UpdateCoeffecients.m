% THIS SCRIPT UPDATES THE NEW COEFFECIENTS :


% USE THE NEW VALUES ON THE APPRXIMATING GRID TO UPDATE THE COEFF:This is
% done using CompEcon's funfitxy routine. This routine internally uses some
% sort of least square fit

    for s_=1:S
            IndxUnSolved_s=find(~(ExitFlag==1));
    IndxSolved_s=(s_-1)*GridSize/S+find(ExitFlag((s_-1)*GridSize/S+1:s_*GridSize/S)==1);
        cNew(s_,:) = funfitxy(V(s_),domain(IndxSolved_s,1:2),VNew(IndxSolved_s)');
        for s=1:S
        Coeff_xhat(s_,s,:) = funfitxy(xhat(s_,s),domain(IndxSolved_s,1:2), PolicyRulesStore(IndxSolved_s,end-S+s));
        Coeff_Rhat(s_,s,:) = funfitxy(Rhat(s_,s),domain(IndxSolved_s,1:2),PolicyRulesStore(IndxSolved_s,end-2*S+s));
        end
   end
        
   % STORE THE DIFF IN COEFF
    cdiff(iter,:)=sum(abs(c-cNew))';
    cOld=c;
    % UPDATE THE COEFF USING A WEIGTED AVERAGE SCHEME
    % update the guess by taking a weighted average of the old and new
    % coeffecients
    c=cNew*Para.grelax+(1-Para.grelax)*cOld;
    
    %ERROR IN SUP NORM
    Error=0;
    for s=1:S
        IndxSolved_s=(s-1)*GridSize/S+find(ExitFlag((s-1)*GridSize/S+1:s*GridSize/S)==1);
    Error=Error+max(abs(VNew(IndxSolved_s)'-funeval(cOld(s,:)',V(s),domain(IndxSolved_s,1:2))));
    end
    ErrorInSupNorm(iter-1)=Error
    disp('Completed Iteration No - ')
    disp(iter)
    toc
    % SAVE THE NEW COEFF
    if mod(iter,1)==0
        save([ Para.datapath Para.StoreFileName] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','domain','Para','V','xhat','Coeff_xhat','Rhat','Coeff_Rhat');
    end
    
