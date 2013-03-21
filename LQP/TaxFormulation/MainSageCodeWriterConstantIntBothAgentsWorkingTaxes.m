    % var('ss_l1,ss_c1,ss_c2,ss_b2,ss_l2,ss_g,ss_Q,g,c1,c2,l1,l2,b2,c1to,c2to,b)
    % 2to,l1to,l2to,sigma,gamma,beta,theta_1,theta_2')
    
    
 clear all
EndoVarList={'tau','T','c1','c2','b2','Q'}
ShocksList={'g','theta_1','theta_2'}

%VarList={'c1','c2'}
ParamList={'sigma','gamma','beta','alpha_1','alpha_2'};
EqName{1}='I'
ExConstr{1}='ss_lambda_I*( ss_g*exp(g)+2*ss_T*exp(T)+ss_Q*exp(Q)*ss_b2*exp(b2)-ss_tau*exp(tau)*(((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_1^(1/gamma)*exp(theta_1/gamma)*ss_c1^(-sigma/gamma)*exp(c1*(-sigma/gamma)))*ss_theta_1*exp(theta_1)+ ((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_2^(1/gamma)*exp(theta_2/gamma)*ss_c2^(-sigma/gamma)*exp(c2*(-sigma/gamma)))*ss_theta_2*exp(theta_2))-ss_b2*exp(b2_))'
EqName{2}='R'
ExConstr{2}='ss_lambda_R*( ss_c1*exp(c1)+ss_c2*exp(c2)+ss_g*exp(g)-  ss_theta_1*exp(theta_1)*((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_1^(1/gamma)*exp(theta_1/gamma)*ss_c1^(-sigma/gamma)*exp(c1*(-sigma/gamma)))- ss_theta_2*exp(theta_2)*((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_2^(1/gamma)*exp(theta_2/gamma)*ss_c2^(-sigma/gamma)*exp(c2*(-sigma/gamma))) )'
EqName{3}='B2'
ExConstr{3}='ss_lambda_B2*(ss_c2*exp(c2)-ss_T*exp(T)-ss_c2^(sigma)*exp(sigma*c2)*((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_2^(1/gamma)*exp(theta_2/gamma)*ss_c2^(-sigma/gamma)*exp(c2*(-sigma/gamma)))^(1+gamma) )'
EqName{4}='E1'
ExConstr{4}='ss_lambda_E1*( beta*ss_c1^(-sigma)*exp(-sigma*c1)-ss_Q*exp(Q_)*ss_c1^(-sigma)*exp(-sigma*c1_)  )'
EqName{5}='E2'
ExConstr{5}='ss_lambda_E2*(  beta*ss_c2^(-sigma)*exp(-sigma*c2)-ss_Q*exp(Q_)*ss_c2^(-sigma)*exp(-sigma*c2_)  )'
EqName{6}='O'
ExConstr{6}='alpha_1*((ss_c1*exp(c1))^(1-sigma)/(1-sigma)-((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_1^(1/gamma)*exp(theta_1/gamma)*ss_c1^(-sigma/gamma)*exp(c1*(-sigma/gamma)))^(1+gamma)/(1+gamma))+alpha_2*((ss_c2*exp(c2))^(1-sigma)/(1-sigma)-((1-ss_tau*exp(tau))^(1/gamma)*ss_theta_2^(1/gamma)*exp(theta_1/gamma)*ss_c2^(-sigma/gamma)*exp(c2*(-sigma/gamma)))^(1+gamma)/(1+gamma)) '
EqName{7}='QQ'
ExConstr{7}='ss_lambda_QQ*(ss_Q*exp(Q)-beta)'


    fid = fopen('SageCodeTaxes.txt', 'wt');


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
    [ SSExpression,EqExpression,LinearApprox,QuadApprox,Dy,Dycheck,Dxi,Dyy,Dycheckycheck,Dyycheck,Dyxi,Dycheckxi,Dxixi] = getSageCode( EqName{i}, ExConstr{i}, EndoVarList,ShocksList)
    
    fprintf(fid,'\n')
    fprintf(fid,EqExpression)
        fprintf(fid,'\n')
    fprintf(fid,SSExpression)

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
