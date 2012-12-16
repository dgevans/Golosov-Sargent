% THIS SCRIPT UPDATES THE NEW COEFFECIENTS :


% USE THE NEW VALUES ON THE APPRXIMATING GRID TO UPDATE THE COEFF:This is
% done using CompEcon's funfitxy routine. This routine internally uses some
% sort of least square fit

    cNew(1,:)=funfitxy(V(1),domain(IndxSolved_1,1:2),VNew(IndxSolved_1)' );
    cNew(2,:)=cNew(1,:);
        
   % STORE THE DIFF IN COEFF
    cdiff(iter,:)=sum(abs(c-cNew))';
    cOld=c;
    % UPDATE THE COEFF USING A WEIGTED AVERAGE SCHEME
    % update the guess by taking a weighted average of the old and new
    % coeffecients
    c=cNew*Para.grelax+(1-Para.grelax)*cOld;
    
    %ERROR IN SUP NORM
    ErrorInSupNorm(iter-1)=max(abs(VNew(IndxSolved_1)'-funeval(cOld(1,:)',V(1),domain(IndxSolved_1,1:2))));
    disp('Completed Iteration No - ')
    disp(iter)
    toc
    % SAVE THE NEW COEFF
    if mod(iter,1)==0
        save([ Para.datapath  'c_' num2str(iter) '.mat' ] , 'c','ErrorInSupNorm','cdiff','IndxSolved','IndxUnSolved','PolicyRulesStore','VNew','domain','Para','V');
    end
    