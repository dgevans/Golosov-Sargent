function [ cNew ] = FitConcaveValueFunctionBEGS(V,VNew,xR)
%FitConcaveValueFunction : This function updates the coeffecients using a
%least square fit imposing concavity with respect to x,R at randomly
%choosen 20% points
%   Detailed explanation goes here
[cguess]=funfitxy(V,xR,VNew);

ShapeTestPoints=xR(1:round(prod(V.n)/10):end,:);
LSRes=@(c) (-sum((funeval(c,V,xR)-VNew).^2))^.5;
NumShapeTestPoints=size(ShapeTestPoints,1);
for ctr=1:NumShapeTestPoints
ShapeConstraintsx(ctr,:)=funbas(V,ShapeTestPoints(ctr,:),[2 0]);
ShapeConstraintsR(ctr,:)=funbas(V,ShapeTestPoints(ctr,:),[0 2]);
end


A=[ShapeConstraintsx;ShapeConstraintsR];
b=zeros(NumShapeTestPoints*2,1);
%A=[];
%b=[];
Aeq=[];
beq=[];
opts = optimset('Algorithm', 'interior-point', 'Display','on');    
cNew = fmincon(@(c) LSRes(c),cguess,A,b,Aeq,beq,[],[],[],opts);




end
