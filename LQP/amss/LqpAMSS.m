    % var('ss_l1,ss_c1,ss_c2,ss_b2,ss_l2,ss_g,beta,g,c1,c2,l1,l2,b2,c1to,c2to,b)
    % 2to,l1to,l2to,sigma,gamma,beta,theta_1,theta_2')
    clear all
    EndoVarList={'c','l','b','Q'}
    ShocksList={'g'}
    %VarList={'c1','c2'}
    ParamList={'sigma','theta','theta_2','gamma','beta'};
    EqName{1}='I'
    ExConstr{1}='ss_lambda_I*(ss_c*exp(c)+ss_b*ss_Q*exp(b+Q)-ss_l^(1+gamma)*ss_c^(sigma)*exp(sigma*c+(1+gamma)*l)-ss_b*exp(b_) )'
    EqName{2}='B'
    ExConstr{2}='ss_lambda_B*(ss_Q*ss_c^(-sigma)*exp(-sigma*c_+Q_)-beta*ss_c^(-sigma)*exp(-sigma*c) )'
    EqName{3}='R'
    ExConstr{3}='ss_lambda_R*(ss_c*exp(c)+ss_g*exp(g)-theta*ss_l*exp(l))'
    EqName{4}='O'
    ExConstr{4}='((ss_c*exp(c))^(1-sigma)/(1-sigma)- (ss_l*exp(l))^(1+gamma)/(1+gamma) )'
    
    fid = fopen('SageCodeAMSSEconomy.txt', 'wt');

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
