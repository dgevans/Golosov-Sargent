function [c1] = getValueC1(x,R,s,Para )
%GETVALUEC1 Summary of this function goes here
% This function solves the Implementability condition for c1. since there
% are two roots it uses the root that gives the highest utility

%   Detailed explanation goes here
    
    n1=Para.n1;
    n2=Para.n2;
    alpha_1=Para.alpha_1;
    alpha_2=Para.alpha_2;
    g=Para.g(s);
    theta_1=Para.theta_1;
    theta_2=Para.theta_2;
    psi=Para.psi;
    sigma=Para.sigma;
    MaxResources=theta_1*n1+theta_2*n2-g;
    %First compute the low root
     %options = optimset('Display','iter');
     options = optimset('Display','off');
     [c1,~,exitflagB,~]= fzero(@(c1) SolveImpCons(c1,R,x,s,Para),(MaxResources./(1+R^(-1/sigma)))/2,options);
%        %SolveImpCons(2,R,x,Para)     
%      % fplot(@(c1) SolveImpCons(c1,R,s,x,Para), [1.6 5])
%      c2=R^(-1)*c1;
%     TotalResources=(c1*n1+c2*n2+g);
% FF=R*theta_2/theta_1;
% DenL2=n1*theta_1*FF+theta_2*n2;
% l2=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
% l1= 1-FF*(1-l2);
% 
% 
%     
%     v1=alpha_1*uBGP(c1,l1,psi)+alpha_2*uBGP(c2,l2,psi);
%     c1_1 = c1;
%     % if the solution is not acceptable put c1=-1 and v1 =-Inf
   if (~(exitflagB==1)|| c1 <0.0001 || isnan(c1))
     c1=-1;

%     v1=-Inf;
     end
% 
%     %high root
%     
%      [c1,~,exitflagB,~]= fzero(@(c1) SolveImpCons(c1,R,x,s,Para),10,options);
%          c2=R^(-1)*c1;
%      TotalResources=(c1*n1+c2*n2+g);
% FF=R*theta_2/theta_1;
% DenL2=n1*theta_1*FF+theta_2*n2;
% l2=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
% l1= 1-FF*(1-l2);
%     v2=alpha_1*uBGP(c1,l1,psi)+alpha_2*uBGP(c2,l2,psi);
%     c1_2 = c1;
%        % if the solution is not acceptable put c1=-1 and v1 =-Inf
%      if (~(exitflagB==1)|| c1 <0.0001||isnan(c1))
%                 c1=-1;
%                 c1_2 = c1;
%                 v2=-Inf;
%        end
%        
%     
%     if( v1 > v2)
%         c1 = c1_1;
%         U = v1;
%     else
%         c1 = c1_2;
%         U = v2;
%     end
% end
% 
