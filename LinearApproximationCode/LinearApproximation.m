function [ A XSS, B, BS ] = LinearApproximation( Para)
%LINEARAPPROXIMATION Function takes a linear approximation to steady state
%   Detailed explanation goes here

    [x,R,PolicyRules] = findSteadyState(0,3,Para);

    XSS = [PolicyRules(1:8) PolicyRules(9) PolicyRules(9) PolicyRules(10) PolicyRules(10)...
       PolicyRules(11) PolicyRules(11) PolicyRules(12:18)];
   
    [VRSS,VxSS] = computeVSS(XSS,Para);
    
    YSS = [R x  VRSS(1) VxSS(1) VRSS(2) VxSS(2)];
    
    f = @(X,Y) FOCResiduals(X,Y(1),Y(2),[Y(3) Y(5)],[Y(4) Y(6)],Para);

    [Dg1 Dg2] = fDerivative(f,XSS,YSS);
    
    DV = computeDV(XSS,Para);
    
    user.Dg1 = Dg1;
    user.Dg2 = Dg2;
    user.DV = DV;
    diff = 1;
    for i = 1:100
        [C, fvec, ~, ifail] = c05qb(@MatrixEquationNAG,randn(42,1),'n',int64(42),'xtol',1e-10,'user',user);
        if(ifail == 0)
            diffnew = norm(fvec,Inf);
            if(diffnew < diff)
                diff = diffnew
                A = reshape(C,21,2);
            end
        end
        
    end
    if(diff == 1)
        throw(MException('Could Not Find Root'));
    end
    
    P = Para.P(1,:);
    B{1} = A([9,11],:);
    B{2} = A([10,12],:);
    fSigma = @(Sigma) P(1)*B{1}*Sigma*B{1}' + P(2)*B{2}*Sigma*B{2}';
    I = eye(4);
    BS = zeros(4);
    for i = 1:4
        BS(:,i) = reshape(fSigma(reshape(I(:,i),2,2)),4,1);
    end
    
   
end

function [fvec, user, iflag] = MatrixEquationNAG(n, x, user, iflag)
    M = reshape(x,21,2);
    Dg1= user.Dg1;
    Dg2 = user.Dg2;
    DV = user.DV;
    
    fvec = reshape(MatrixEquation(Dg1,Dg2,DV,M)-M,42,1);
    
end

function [VRSS,VxSS] = computeVSS(XSS,Para)
    sigma = Para.sigma;
    beta = Para.beta;
    psi = Para.psi;
    P = Para.P(1,:);
   
    c1 = XSS(1:2);
    c2 = XSS(3:4);
    lambda = XSS(15);
    mu = XSS(13:14);

    uc1 = psi*c1.^(-sigma);
    uc2 = psi*c2.^(-sigma);

    VxSS = dot(uc2,mu.*P)/(beta*dot(uc2,P))*ones(1,2);
    VRSS = -lambda*dot(uc1,P)*ones(1,2);
    

end


function [Dg1 Dg2] = fDerivative(f,XSS,YSS)
    IX = eye(21);
    IY = eye(6);
    
    f0 = f(XSS,YSS)';
    DX = zeros(21);
    h = 1e-6;
    for i = 1:21
        X = XSS +h*IX(i,:);
        DX(:,i) = (f(X,YSS)'-f0)/h;
    end
    DY = zeros(21,6);
    for i = 1:6
        Y = YSS +h*IY(i,:);
        DY(:,i) = (f(XSS,Y)'-f0)/h;
    end
    
    DG = -DX\DY;
    
    Dg1 = DG(:,1:2);
    Dg2 = DG(:,3:6);
end

function [Aprime] = MatrixEquation(Dg1,Dg2,DV,A)

    IX = eye(21);
    eRx = [IX(9,:); IX(11,:); IX(10,:); IX(12,:)];
    
    Aprime = Dg1 + Dg2*[DV zeros(2,21);zeros(2,21) DV]*[A zeros(21,2);zeros(21,2) A]*eRx*A;
end

function [DV] = computeDV(XSS,Para)
    beta = Para.beta;
    P = Para.P(1,:);
    U = Para.U;
    
    c1= XSS(1:2);
    c2 = XSS(3:4);
    mu = XSS(13:14);
    lambda = XSS(15);
    [~,uc1,~,ucc1] = U(c1,0.5*ones(1,2));
    [~,uc2,~,ucc2] = U(c2,0.5*ones(1,2));
    Euc1 = dot(P,uc1);
    Euc2 = dot(P,uc2);
    Euc2_mu = dot(P,uc2.*mu);

    
    DV = zeros(2,21);
    
    DVR_dc1 = -lambda*P.*ucc1;
    DVR_dlambda = -Euc1;
    DV(1,1:2) = DVR_dc1;
    DV(1,15) = DVR_dlambda;
    
    DVx_dc2 = mu.*P.*ucc2/(beta*Euc2) - (P.*ucc2*Euc2_mu)/(beta*Euc2^2);
    DVx_dmu = P.*uc2/(beta*Euc2);
    
    DV(2,3:4) = DVx_dc2;
    DV(2,13:14) = DVx_dmu;
end
