function [res] = SolveNoShockProblem(xState,Para)
u2btild_=xState(1);
R_=xState(2);
s_=xState(3);
                    c1_=max(getValueC1(u2btild_,R_,s_,Para ),.0001);
                    
                  
                    if c1_<.001
                        ExitFlagT=0;
                    else
                        ExitFlagT=1;
                    end
                    % compute c2
                    c2_=R_^(-1)*c1_;
                    TotalResources=(c1_*Para.n1+c2_*Para.n2+Para.g(s_));
                    FF=R_*Para.theta_2/Para.theta_1;
                    DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
                    l2_=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
                    if Para.theta_2==0
    l1_=TotalResources/Para.theta_1;
l2_=0;
else
l1_= 1-FF*(1-l2_);
end

                    u2btildPrime_=u2btild_;
                    Vobj=(Para.alpha_1*uBGP(c1_,l1_,Para.psi)+Para.alpha_2*uBGP(c2_,l2_,Para.psi))/(1-Para.beta);
PolicyRules=[c1_ c1_ c2_];
res.PolicyRules=PolicyRules;
res.Value=Vobj;
res.exitflag=ExitFlagT;
end

