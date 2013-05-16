function [ domain, c, PolicyRulesStore] = InitializeCoeffWithCompleteMarkets( Para, V)
S = length(Para.P);
xGrid=Para.xGrid;
RGrid=Para.RGrid;
xGridSize=Para.xGridSize;
RGridSize=Para.RGridSize;
flagIID =min(Para.P(:,1))==max(Para.P(:,1));
PolicyRulesStore=[];
domain=[];
VNew0=[];
          
for s_=1:S
        
        for xctr=1:Para.xGridSize
            for Rctr=1:Para.RGridSize
                
                x=xGrid(xctr);
                R=RGrid(Rctr);
                if Rctr+xctr+s_<6
                      cRat = R^(-1/Para.sigma);
                c1 = (0.8*(Para.n1*Para.theta_1+Para.n2*Para.theta_2)-Para.g)./(Para.n1+cRat*Para.n2);
                c1 = c1(:)';
                c2_ = cRat*c1; c2_(S) = []; 
                options = optimset('Display','off');
                warning off
                [xSS,~,exitFlag] = fsolve(@(z) SteadyStateResiduals(z,x,R,Para,s_),[c1 c2_],options);
                [res, c1, c2, l1, l2] = SteadyStateResiduals(xSS,x,R,Para,s_);
                z0=[c1 c2 l1 l2 zeros(1,3*S) 0];

                end
                
                
                [ PolicyRule,VNew,z0,ifail ] = getCompleteMarketSolution(x,R,s_,Para,z0);
                
                domain=[domain;x R s_];
                PolicyRulesStore=[PolicyRulesStore;PolicyRule];
                VNew0=[VNew0;VNew];
                if ~(ifail==0)
                    disp([x R s_])
                end
                    
                %end
            end
        end
        c(s_,:)=funfitxy(V(s_),domain((s_-1)*xGridSize*RGridSize+1:s_*xGridSize*RGridSize,1:2),VNew0((s_-1)*xGridSize*RGridSize+1:s_*xGridSize*RGridSize ));
        if  flagIID
    
        if ~(s_==1)
            for s=2:S
                c(s,:)=c(1,:);
                domain=[domain; domain((s_-1)*xGridSize*RGridSize+1:s_*xGridSize*RGridSize,:)+repmat([1 0 0],xGridSize*RGridSize,1)];
                PolicyRulesStore=[PolicyRulesStore;PolicyRulesStore((s_-1)*xGridSize*RGridSize+1:s_*xGridSize*RGridSize,:)];
               end
            
        break;
        end
        
        end
end
        
   
end

