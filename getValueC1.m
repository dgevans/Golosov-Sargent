function [c1] = getValueC1(u2btild,R,s,Para )
%GETVALUEC1 Summary of this function goes here
% This function solves the Implementability condition for c1. since there
% are two roots it uses the root that gives the highest utility

%   Detailed explanation goes here
    
    n1=Para.n1;
    n2=Para.n2;
    alpha_1=Para.alpha_1;
    alpha_2=Para.alpha_2;
    g=Para.g;
    theta_1=Para.theta_1(s);
    theta_2=Para.theta_2(s);
    psi=Para.psi;
    
    
    FF=R*theta_2/theta_1;
cUpperBound=(theta_1*n1*FF+theta_2*n2-theta_1*n1*(FF-1)-g)/(n1+n2/R);

    %First compute the low root
     options = optimset('Display','off');
     [c1,~,exitflagB,~]= fzero(@(c1) SolveImpCons(c1,R,u2btild,s,Para),cUpperBound/2,options);
     % if the solution is not acceptable put c1=-1 and v1 =-Inf
   if (~(exitflagB==1)|| c1 <0.0001 || isnan(c1))
     c1=-1;

     end
