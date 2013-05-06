function [ cNew ] = SmoothPasting(V,VNew,x,DerivativeMatch,DerivativeVals,ShapeTestPoints)
%   Detailed explanation goes here
[cguess]=funfitxy(V,x,VNew);
LSRes=@(c) (sum((funeval(c,V,x)-VNew).^2))^.5;
NumShapeTestPoints=length(ShapeTestPoints);
MontonicityConstraint=[];
for ctr=1:NumShapeTestPoints
mc=funbasx(V,ShapeTestPoints(ctr,:),[1],'expanded');
MontonicityConstraint(ctr,:)=mc.vals{1};
end
DerivativeMatch_c=funbasx(V,DerivativeMatch,[1],'expanded');

SmoothPastingConstraint=DerivativeMatch_c.vals{1};

A=[MontonicityConstraint];
b=[zeros(length(MontonicityConstraint),1)];
%A=[];
%b=[];
Aeq=[SmoothPastingConstraint];
beq=[DerivativeVals];
opts = optimset('Display','off','Algorithm','interior-point');    
[cNew ,fval,exitflag] = fmincon(@(c) LSRes(c),cguess,A,b,Aeq,beq,[],[],[],opts);

cNew=cNew';


end
