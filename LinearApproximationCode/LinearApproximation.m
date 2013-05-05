function [ A, XSS, B, BS , DX,Dy1,Dy2,DV] = LinearApproximation( Para)
%LINEARAPPROXIMATION Function takes a linear approximation to steady state
%   Detailed explanation goes here

    [x,R,PolicyRules] = findSteadyState(0,3,Para);
if min(Para.theta_2)==0
    Para.dim_foc=17;
    Para.offset_l2=2;
    
    XSS = [PolicyRules(1:8-2) PolicyRules(9-2) PolicyRules(9-2) PolicyRules(10-2) PolicyRules(10-2)...
       PolicyRules(11-2) PolicyRules(11-2) PolicyRules(12-2:18-4)];
   
else
    % if length(PolicyRules)==14
%     PolicyRules(1:6)=PolicyRules(1:6);
%     PolicyRules(7:8)=[0 0];
%     PolicyRules(9:12)=PolicyRules(7:10);
%     PolicyRules(13:14)=[0 0];
%     PolicyRules(15:18)=PolicyRules(11:14);
% end
    Para.dim_foc=21;
    
Para.offset_l2=0;    
    XSS = [PolicyRules(1:8) PolicyRules(9) PolicyRules(9) PolicyRules(10) PolicyRules(10)...
       PolicyRules(11) PolicyRules(11) PolicyRules(12:18)];
end
    [VxSS,VRSS] = computeVSS(XSS,Para);
    
    YSS = [x R  VxSS(1) VRSS(1) VxSS(2) VRSS(2)];
    
    f = @(X,Y,Para) FOCResiduals(X,Y(1),Y(2),[Y(3) Y(5)],[Y(4) Y(6)],Para);

    [DX Dy1 Dy2] = fDerivative(f,XSS,YSS,Para);
    
    DV = computeDV(XSS,Para);
    
    user.Dy1 = Dy1;
    user.Dy2 = Dy2;
    user.DV = DV;
    user.DX = DX;
    user.Para=Para;
    diff = 1;
    for i = 1:100
        [C, fvec, ~, ifail] = c05qb(@MatrixEquationNAG,randn(Para.dim_foc*2,1),'n',int64(Para.dim_foc*2),'xtol',1e-10,'user',user);
        if(ifail == 0)
            diffnew = norm(fvec,Inf);
            if(diffnew < diff)
                diff = diffnew
                A = reshape(C,Para.dim_foc,2);
            end
        end
        
    end
    if(diff == 1)
        throw(MException('Could Not Find Root'));
    end
    
    P = Para.P(1,:);
    B{1} = A([9-Para.offset_l2,11-Para.offset_l2],:);
    B{2} = A([10-Para.offset_l2,12-Para.offset_l2],:);
    fSigma = @(Sigma) P(1)*B{1}*Sigma*B{1}' + P(2)*B{2}*Sigma*B{2}';
    I = eye(4);
    BS = zeros(4);
    for i = 1:4
        BS(:,i) = reshape(fSigma(reshape(I(:,i),2,2)),4,1);
    end
    
   
end

function [fvec, user, iflag] = MatrixEquationNAG(n, x, user, iflag)
    M = reshape(x,user.Para.dim_foc,2);
    Dy1= user.Dy1;
    Dy2 = user.Dy2;
    DX = user.DX;
    DV = user.DV;
    
    fvec = reshape(MatrixEquation(DX,Dy1,Dy2,DV,M,user.Para),user.Para.dim_foc*2,1);
    
end


function [VxSS,VRSS] = computeVSS(XSS,Para)
    beta = Para.beta;
    P = Para.P(1,:);
    U = Para.U;
   
    c1 = XSS(1:2);
    c2 = XSS(3:4);
    lambda = XSS(15-Para.offset_l2);
    mu = XSS(13-Para.offset_l2:14-Para.offset_l2);

    
%    [~,uc1,] = U(c1,0.5*ones(1,2),Para);
 %   [~,uc2,] = U(c2,0.5*ones(1,2),Para);
    [~,uc1,] = U(c1,0.5*ones(1,2));
    [~,uc2,] = U(c2,0.5*ones(1,2));

    
    VxSS = dot(uc2,mu.*P)/(beta*dot(uc2,P))*ones(1,2);
    VRSS = -lambda*dot(uc1,P)*ones(1,2);
    

end


function [DX,Dy1 Dy2] = fDerivative(f,XSS,YSS,Para)
    IX = eye(Para.dim_foc);
    IY = eye(6);
    
    f0 = f(XSS,YSS,Para)';
    DX = zeros(Para.dim_foc);
    h = 1e-6;
    for i = 1:Para.dim_foc
        X = XSS +h*IX(i,:);
        DX(:,i) = (f(X,YSS,Para)'-f0)/h;
    end
    DY = zeros(Para.dim_foc,6);
    for i = 1:6
        Y = YSS +h*IY(i,:);
        DY(:,i) = (f(XSS,Y,Para)'-f0)/h;
    end
   
    
    Dy1 = DY(:,1:2);
    Dy2 = DY(:,3:6);
end

function [ret] = MatrixEquation(DX,Dy1,Dy2,DV,A,Para)

    IX = eye(Para.dim_foc);
    exR = [IX(9- Para.offset_l2,:); IX(11- Para.offset_l2,:); IX(10- Para.offset_l2,:); IX(12- Para.offset_l2,:)];
    
    ret = DX*A +Dy1 + Dy2*[DV zeros(2,Para.dim_foc);zeros(2,Para.dim_foc) DV]*[A zeros(Para.dim_foc,2);zeros(Para.dim_foc,2) A]*exR*A;
end

function [DV] = computeDV(XSS,Para)
    beta = Para.beta;
    P = Para.P(1,:);
    U = Para.U;
    
    c1= XSS(1:2);
    c2 = XSS(3:4);
    mu = XSS(13-Para.offset_l2:14-Para.offset_l2);
    lambda = XSS(15-Para.offset_l2);
%    [~,uc1,] = U(c1,0.5*ones(1,2),Para);
 %   [~,uc2,] = U(c2,0.5*ones(1,2),Para);
[~,uc1,~,ucc1,~] = U(c1,0.5*ones(1,2));
    [~,uc2,~,ucc2,~] = U(c2,0.5*ones(1,2));
   
    Euc1 = dot(P,uc1);
    Euc2 = dot(P,uc2);
    Euc2_mu = dot(P,uc2.*mu);

    
    DV = zeros(2,Para.dim_foc);
        
    DVx_dc2 = mu.*P.*ucc2/(beta*Euc2) - (P.*ucc2*Euc2_mu)/(beta*Euc2^2);
    DVx_dmu = P.*uc2/(beta*Euc2);
    
    DV(1,3:4) = DVx_dc2;
    DV(1,13-Para.offset_l2:14-Para.offset_l2) = DVx_dmu;
    
    
    DVR_dc1 = -lambda*P.*ucc1;
    DVR_dlambda = -Euc1;
    DV(2,1:2) = DVR_dc1;
    DV(2,15-Para.offset_l2) = DVR_dlambda;

end
