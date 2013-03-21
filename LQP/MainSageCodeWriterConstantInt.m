% var('ss_l1,ss_c1,ss_c2,ss_b2,ss_l2,ss_g,ss_Q,g,c1,c2,l1,l2,b2,c1to,c2to,b)
% 2to,l1to,l2to,sigma,gamma,beta,theta_1,theta_2')
clear all
EndoVarList={'tau','T','l1','l2','b'}
ShocksList={'g','theta_1','theta_2'}
%VarList={'c1','c2'}
ParamList={'sigma','gamma','beta','alpha_1','alpha_2'};

EqName{1}='GBC'
ExConstr{1}='ss_lambda_GBC*(ss_g*exp(g)+ss_T*exp(T)+beta*ss_b*exp(b)-ss_tau*ss_y*exp(tau+y)-ss_b*exp(b_))'
EqName{2}='R'
ExConstr{2}='ss_lambda_R*(theta_1^(1+gamma)/(sigma) *(1-ss_tau*exp(tau))^(1/sigma)*(ss_y*exp(y))^(-gamma/sigma) + ss_T*exp(T)+ ss_g*exp(g)-ss_y*exp(y))'
EqName{3}='O'
ExConstr{3}='alpha_1*((theta_1^(1+gamma)/(sigma) *(1-ss_tau*exp(tau))^(1/sigma)*(ss_y*exp(y))^(-gamma/sigma) )^(1-sigma)/(1-sigma)- (ss_y*exp(y)/theta_1)^(1+gamma)/(1+gamma) )+ alpha_2*((ss_T*exp(T))^(1-sigma)/(1-sigma))'


fid = fopen('SageCodeConstInt.txt', 'wt');


n=length(EndoVarList);
ns=length(ShocksList);

for i=1:length(ParamList)
fprintf(fid,'\n')
fprintf(fid,['var ('''   ParamList{i} ''')'])

end


for i=1:n
VarList{i}=[EndoVarList{i}];
fprintf(fid,'\n')
fprintf(fid,['var ('''  'ss_' VarList{i} ''')'])

end
for i=1:ns
fprintf(fid,'\n')
fprintf(fid,['var ('''  'ss_' ShocksList{i} ''')'])

end



for i=1:n
VarList{i}=[EndoVarList{i} '_'];
fprintf(fid,'\n')
fprintf(fid,['var ('''  VarList{i} ''')'])
end
for i=n+1:2*n
VarList{i}=[EndoVarList{i-n}];
fprintf(fid,'\n')
fprintf(fid,['var ('''  VarList{i} ''')'])

end

for i=2*n+1:2*n+ns
VarList{i}=[ShocksList{i-2*n}];
fprintf(fid,'\n')
fprintf(fid,['var ('''  VarList{i} ''')'])
end


for i=1:length(EqName)

fprintf(fid,'\n')

fprintf(fid,['var ('''  'ss_lambda_' EqName{i} ''')'])
[ EqExpression,LinearApprox,QuadApprox,Dy,Dycheck,Dxi,Dyy,Dycheckycheck,Dyycheck,Dyxi,Dycheckxi,Dxixi] = getSageCode( EqName{i}, ExConstr{i}, EndoVarList,ShocksList)
fprintf(fid,'\n')
fprintf(fid,EqExpression)
fprintf(fid,'\n')
fprintf(fid,LinearApprox)
fprintf(fid,'\n')
fprintf(fid,QuadApprox)
fprintf(fid,'\n')

fprintf(fid,['# Linear Terms']);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dy = ' Dy]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dycheck = ' Dycheck]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dxi = ' Dxi]);

fprintf(fid,'\n')
fprintf(fid,['# Quadratic Terms']);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dyy = ' Dyy]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dycheckycheck = ' Dycheckycheck]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dyycheck = ' Dyycheck]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dyxi = ' Dyxi]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dycheckxi= ' Dycheckxi]);
fprintf(fid,'\n')
fprintf(fid,[EqName{i} 'Dxixi = ' Dxixi]);
fprintf(fid,'\n')
fprintf(fid,'\n')


%,LinearApprox,,Dy,Dycheck,Dxi,Dyy,Dycheckycheck,Dyycheck,Dyxi,Dycheckxi,Dxixi)

end
fclose(fid)
