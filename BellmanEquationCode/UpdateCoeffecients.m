% THIS SCRIPT UPDATES THE NEW COEFFECIENTS :


% USE THE NEW VALUES ON THE APPRXIMATING GRID TO UPDATE THE COEFF:This is
% done using CompEcon's funfitxy routine. This routine internally uses some
% sort of least square fit

    cNew(1,:)=funfitxy(V(1),domain(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
        Coeff_xhat(1,:)=funfitxy(xhat(1),domain(IndxSolved_1,1:2), PolicyRulesStore(IndxSolved_1,end-S+1));
        Coeff_Rhat(1,:)=funfitxy(Rhat(1),domain(IndxSolved_1,1:2), PolicyRulesStore(IndxSolved_1,end-2*S+1));
    for s=2:S
        cNew(s,:) = cNew(1,:);
        Coeff_xhat(s,:)=funfitxy(xhat(s),domain(IndxSolved_1,1:2), PolicyRulesStore(IndxSolved_1,end-S+s));
        Coeff_Rhat(s,:)=funfitxy(Rhat(s),domain(IndxSolved_1,1:2), PolicyRulesStore(IndxSolved_1,end-2*S+s));
    end
        
   % STORE THE DIFF IN COEFF
    cdiff(iter,:)=sum(abs(c-cNew))';
    cOld=c;
    % UPDATE THE COEFF USING A WEIGTED AVERAGE SCHEME
    % update the guess by taking a weighted average of the old and new
    % coeffecients
    c=cNew*Para.grelax+(1-Para.grelax)*cOld;
    
    %ERROR IN SUP NORM
    ErrorInSupNorm(iter-1)=max(abs(VNew(IndxSolved_1)'-funeval(cOld(1,:)',V(1),domain(IndxSolved_1,1:2))))
    disp('Completed Iteration No - ')
    disp(iter)
    toc
    % SAVE THE NEW COEFF
    if mod(iter,5)==0
        save([ Para.datapath Para.StoreFileName] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','domain','Para','V','xhat','Coeff_xhat','Rhat','Coeff_Rhat');
    end
    
