function [ cNew ] = SmoothPasting(V,VNew,x,x_cutoff,ShapeTestPoints)
%   Detailed explanation goes here
[cguess]=funfitxy(V,x,VNew);
LSRes=@(c) (sum((funeval(c,V,x)-VNew).^2))^.5;
NumShapeTestPoints=length(ShapeTestPoints);
MontonicityConstraint=[];
for ctr=1:NumShapeTestPoints
mc=funbasx(V,ShapeTestPoints(ctr,:),[1],'expanded');
MontonicityConstraint(ctr,:)=mc.vals{1};
end
cutoff_der_c=funbasx(V,x_cutoff,[1],'expanded');

SmoothPastingConstraint=cutoff_der_c.vals{1};

A=[MontonicityConstraint];
b=[zeros(length(MontonicityConstraint),1)];
%A=[];
%b=[];
Aeq=[SmoothPastingConstraint];
beq=[0];
opts = optimset('Display','off','Algorithm','interior-point');    
[cNew ,fval,exitflag] = fmincon(@(c) LSRes(c),cguess,A,b,Aeq,beq,[],[],[],opts);

cNew=cNew';


end
