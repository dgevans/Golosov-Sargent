% var('ss_l1,ss_c1,ss_c2,ss_b2,ss_l2,ss_g,ss_Q,g,c1,c2,l1,l2,b2,c1to,c2to,b
% 2to,l1to,l2to,sigma,gamma,beta,theta_1,theta_2')
clear all
EndoVarList={'c1','c2','l1','l2','b2','Q'}
ShocksList={'g'}
%VarList={'c1','c2'}

EqName='1B'
ExConstr='ss_Q*(ss_c1*exp(c1_))^(-sigma)*exp(Q)-beta*(ss_c1*exp(c1))^(-sigma)'



[ EqExpression,LinearApprox,QuadApprox,Dy,Dycheck,Dxi,Dyy,Dycheckycheck,Dyycheck,Dyxi,Dycheckxi,Dxixi] = getSageCode( EqName, ExConstr, EndoVarList,ShocksList)

n=6
ns=1;
