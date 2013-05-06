function [ SSExpression,EqExpression,LinearApprox,QuadApprox,Dy,Dycheck,Dxi,Dyy,Dycheckycheck,Dyycheck,Dyxi,Dycheckxi,Dxixi] = getSageCode( EqName, ExConstr, EndoVarList,ShocksList)
EqExpression=[EqName '=' ExConstr];




SSExpression=['SS_' EqName '=' EqName '(ss_lambda_' EqName '=1,'];
n=length(EndoVarList);
ns=length(ShocksList);
ArgTaylor=[','];
for i=1:n
VarList{i}=[EndoVarList{i} '_'];
ArgTaylor=[ArgTaylor '(' EndoVarList{i} '_,0),'];
SSExpression=[SSExpression EndoVarList{i} '_=0,'];
end
for i=n+1:2*n
VarList{i}=[EndoVarList{i-n}];
ArgTaylor=[ArgTaylor '(' EndoVarList{i-n} ',0),'];
SSExpression=[SSExpression EndoVarList{i-n} '=0,'];
end

for i=2*n+1:2*n+ns
VarList{i}=[ShocksList{i-2*n}];
ArgTaylor=[ArgTaylor '(' ShocksList{i-2*n} ',0),'];
SSExpression=[SSExpression ShocksList{i-2*n} '=0,'];
end
SSExpression(end)=[];
SSExpression=[SSExpression ')'];
LinearApprox= [EqName 'LinearApprox= simplify((taylor(' EqName ArgTaylor '1)))'];
QuadApprox= [EqName 'QuadApprox= simplify((taylor(' EqName ArgTaylor '2)))' '-' EqName 'LinearApprox'];


% Compute the Dy,Dycheck Dxi
Cstr ='vector(['  ; 

for i=1:n
    

    Cstr=[ Cstr ' ' EqName 'LinearApprox.coefficient(' VarList{i} ')'];
    if i<n
        Cstr=[ Cstr ' ,' ];
    end
end

Cstr=[ Cstr '])'];

Dycheck=Cstr;

% compute Dy
Cstr ='vector(['   ;

for i=n+1:2*n
    

    Cstr=[ Cstr ' ' EqName 'LinearApprox.coefficient(' VarList{i} ')'];
    if i<2*n
        Cstr=[ Cstr ' ,' ];
    end
end

Cstr=[ Cstr '])'];

Dy=Cstr;


% compute Dxi
Cstr ='vector(['   ;

for i=2*n+1:2*n+ns
    

    Cstr=[ Cstr ' ' EqName 'LinearApprox.coefficient(' VarList{i} ')'];
    if i<2*n+ns
        Cstr=[ Cstr ' ,' ];
    end
end

Cstr=[ Cstr '])'];

Dxi=Cstr;




% Compute Dyy Dyycheck Dycheckycheck Dyxi Dycheckxi Dxixi

% 1. Dycheckycheck
Cstr ='matrix(['   ;

for i=1:n
    
Cstr=[ Cstr '['    ];
for j=1:n
    
    if ~(i==j)
    Cstr=[ Cstr ' (' EqName 'QuadApprox.coefficient(' VarList{i} ')/2).coefficient(' VarList{j} ')'];
    else
    Cstr=[ Cstr ' ' EqName 'QuadApprox.coefficient(' VarList{i} '^2)'];
    end
    if j<n
    Cstr=[Cstr ',']    ;
    end
end

Cstr=[ Cstr ']'];
if (i<n)
    Cstr=[ Cstr ', '];
else
    Cstr=[ Cstr ' ]) '];
end
end

Dycheckycheck=Cstr;

% 2. compute Dyy

Cstr ='matrix(['   ;

for i=n+1:2*n
    
Cstr=[ Cstr '['    ];
for j=n+1:2*n
    
    if ~(i==j)
    Cstr=[ Cstr ' (' EqName 'QuadApprox.coefficient(' VarList{i} ')/2).coefficient(' VarList{j} ')'];
    else
    Cstr=[ Cstr ' ' EqName 'QuadApprox.coefficient(' VarList{i} '^2)'];
    end
    if j<2*n
    Cstr=[Cstr ',']    ;
    end
end

Cstr=[ Cstr ']'];
if (i<2*n)
    Cstr=[ Cstr ', '];
else
    Cstr=[ Cstr ' ]) '];
end
end

Dyy=Cstr;


% compute Dyycheck
Cstr ='matrix(['   ;

for i=n+1:2*n
    
Cstr=[ Cstr '['    ];
for j=1:n
    
    if ~(i==j)
    Cstr=[ Cstr ' (' EqName 'QuadApprox.coefficient(' VarList{i} ')/2).coefficient(' VarList{j} ')'];
    else
    Cstr=[ Cstr ' ' EqName 'QuadApprox.coefficient(' VarList{i} '^2)'];
    end
    if j<n
    Cstr=[Cstr ',']    ;
    end
end

Cstr=[ Cstr ']'];
if (i<2*n)
    Cstr=[ Cstr ', '];
else
    Cstr=[ Cstr ' ]) '];
end
end

Dyycheck=Cstr;


% compute Dyxi
Cstr ='matrix(['   ;

for i=n+1:2*n
    
Cstr=[ Cstr '['    ];
for j=2*n+1:2*n+ns
    
    if ~(i==j)
    Cstr=[ Cstr ' (' EqName 'QuadApprox.coefficient(' VarList{i} ')/2).coefficient(' VarList{j} ')']
    else
    Cstr=[ Cstr ' ' EqName 'QuadApprox.coefficient(' VarList{i} '^2)']
    end
    if j<2*n+ns
    Cstr=[Cstr ',']    
    end
end

Cstr=[ Cstr ']']
if (i<2*n)
    Cstr=[ Cstr ', '];
else
    Cstr=[ Cstr ' ]) '];
end
end

Dyxi=Cstr;

Cstr ='matrix(['   ;

for i=1:n
    
Cstr=[ Cstr '['    ];
for j=2*n+1:2*n+ns
    
    if ~(i==j)
    Cstr=[ Cstr ' (' EqName 'QuadApprox.coefficient(' VarList{i} ')/2).coefficient(' VarList{j} ')']
    else
    Cstr=[ Cstr ' ' EqName 'QuadApprox.coefficient(' VarList{i} '^2)'];
    end
    if j<2*n+ns
    Cstr=[Cstr ',']    ;
    end
end

Cstr=[ Cstr ']'];
if (i<n)
    Cstr=[ Cstr ', '];
else
    Cstr=[ Cstr ' ]) '];
end
end

Dycheckxi=Cstr;



Cstr ='matrix(['   ;

for i=2*n+1:2*n+ns
    
Cstr=[ Cstr '['    ];
for j=2*n+1:2*n+ns
    
    if ~(i==j)
    Cstr=[ Cstr ' (' EqName 'QuadApprox.coefficient(' VarList{i} ')/2).coefficient(' VarList{j} ')'];
    else
    Cstr=[ Cstr ' ' EqName 'QuadApprox.coefficient(' VarList{i} '^2)'];
    end
    if j<2*n+ns
    Cstr=[Cstr ',']    ;
    end
end

Cstr=[ Cstr ']'];
if (i<2*n+ns)
    Cstr=[ Cstr ', '];
else
    Cstr=[ Cstr ' ]) '];
end
end

Dxixi=Cstr;

end

