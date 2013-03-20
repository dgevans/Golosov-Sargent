function [c1] = getValueC1(x,R,s,Para )
%GETVALUEC1 Summary of this function goes here
% This function solves the Implementability condition for c1.

%RETRIVE THE PARAM    
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
    
    % FIND THE ROOT IN C1(X,R) USING IMPLMENTABILITY CONDITION 
    options = optimset('Display','off');
     [c1,~,exitflagB,~]= fzero(@(c1) SolveImpCons(c1,R,x,s,Para),(MaxResources./(1+R^(-1/sigma)))/2,options);
%     % if the solution is not acceptable put c1=-1 and v1 =-Inf
   if (~(exitflagB==1)|| c1 <0.0001 || isnan(c1))
     c1=-1;

     end
