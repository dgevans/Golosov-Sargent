function [ A XSS B BS Bbar ] = LinearApproximationN( Para)
%LINEARAPPROXIMATION Function takes a linear approximation to steady state
%   Detailed explanation goes here
    N = Para.N;
    S=2;
    [xSS,RSS,PolicyRules] = findSteadyStateN(zeros(N-1,1),repmat(3,N-1,1),Para);
    xSS = repmat(xSS,S,1);
    RSS = repmat(RSS,S,1);
    muSS = repmat(PolicyRules((2*S+2)*N-1:(2*S+3)*N-3)',S,1);
    
    
    XSS = [PolicyRules(1:4*N); reshape(xSS,S*(N-1),1); reshape(RSS,S*(N-1),1);...
       reshape(muSS,S*(N-1),1); PolicyRules((2*S+3)*N-2:end)];
   
   
    [VRSS,VxSS] = computeVSS(XSS,Para);
    
    y1SS = [xSS(1,:) , RSS(1,:)]';
    y2SS = reshape([VxSS';VRSS'],2*S*(N-1),1);
    
    f = @(X,y1,y2) F(X,y1,y2,Para);

    [DX Dy1 Dy2] = fDerivative(f,XSS,y1SS,y2SS);
    
    DV = computeDV(XSS,Para);
    
    user.Dy1 = Dy1;
    user.Dy2 = Dy2;
    user.DV = DV;
    user.DX = DX;
    user.N = N;
    nstate = 2*(N-1);
    diff = 1;
    for i = 1:100
        [C, fvec, ~, ifail] = c05qb(@MatrixEquationNAG,rand(length(XSS)*nstate,1),'n',int64(length(XSS)*nstate),'xtol',1e-10,'user',user);
        if(ifail == 0)
            diffnew = norm(fvec,Inf);
            if(diffnew < diff)
                diff = diffnew
                A = reshape(C,length(XSS),nstate);
            end
        end
        
    end
    if(diff == 1)
        throw(MException('Could Not Find Root'));
    end
    
    P = Para.P(1,:);
    B{1} = A(4*N+1:2:4*N+4*(N-1),:);
    B{2} = A(4*N+2:2:4*N+4*(N-1),:);
    Bbar = P(1)*B{1}+P(2)*B{2};
    nB = length(B{1});
    fSigma = @(Sigma) P(1)*B{1}*Sigma*B{1}' + P(2)*B{2}*Sigma*B{2}';
    I = eye(nB^2);
    BS = zeros(nB^2);
    for i = 1:nB^2
        BS(:,i) = reshape(fSigma(reshape(I(:,i),nB,nB)),nB^2,1);
    end
    
end



function res = F(X,y1,y2,Para)
    S = 2;
    N = Para.N;
    y2 = reshape(y2,2*(N-1),S)';
    Vx = y2(:,1:N-1);
    VR = y2(:,N:2*N-2);
    res = FOCResidualsN(X,y1(1:N-1),y1(N:2*N-2),Vx,VR,Para);
end


function [fvec, user, iflag] = MatrixEquationNAG(n, x, user, iflag)
    nstate = 2*(user.N-1);
    M = reshape(x,n/nstate,nstate);
    Dy1= user.Dy1;
    Dy2 = user.Dy2;
    DV = user.DV;
    DX = user.DX;
    
    fvec = reshape(MatrixEquation(DX,Dy1,Dy2,DV,M,user.N),n,1);
    
end

function [VRSS,VxSS] = computeVSS(X,Para)
    S=2;
    N = Para.N;
    beta = Para.beta;
    P = Para.P(1,:);
    U = Para.U;
   
    ci = reshape(X(1:(N-1)*S),S,N-1);
    X(1:(N-1)*S) = [];
    cN = repmat(reshape(X(1:S),S,1),1,N-1);
    X(1:S) = [];
    X(1:(N-1)*S) = [];
    X(1:S) = [];
    X(1:S*(N-1)) = [];
    X(1:S*(N-1)) = [];
    mu = reshape(X(1:(N-1)*S),S,N-1);
    X(1:S*(N-1)) = [];
    lambda = repmat(reshape(X(1:N-1),1,N-1),S,1);
    
    
    [~,uci] = U(ci,0.5*ones(2,N-1));
    [~,ucN] = U(cN,0.5*ones(2,N-1));
    Emu_ucN = repmat(P * ( mu.*ucN ),S,1);
    Euci = repmat(P * uci,S,1);
    EucN = repmat(P * ucN,S,1);
    
    VxSS = Emu_ucN./(beta*EucN);
    VRSS = -lambda.*Euci;
    

end


function [DX Dy1 Dy2] = fDerivative(f,XSS,y1SS,y2SS)
    nX = length(XSS);
    ny1 = length(y1SS);
    ny2 = length(y2SS);
    IX = eye(nX);
    Iy1 = eye(ny1);
    Iy2 = eye(ny2);
    
    f0 = f(XSS,y1SS,y2SS);
    DX = zeros(nX);
    h = 1e-6;
    for i = 1:nX
        X = XSS +h*IX(:,i);
        DX(:,i) = (f(X,y1SS,y2SS)-f0)/h;
    end
    Dy1 = zeros(nX,ny1);
    for i = 1:ny1
        y1 = y1SS +h*Iy1(:,i);
        Dy1(:,i) = (f(XSS,y1,y2SS)-f0)/h;
    end
    Dy2 = zeros(nX,ny2);
    for i = 1:ny2
        y2 = y2SS +h*Iy2(:,i);
        Dy2(:,i) = (f(XSS,y1SS,y2)-f0)/h;
    end
end

function [res] = MatrixEquation(DX,Dy1,Dy2,DV,A,N)
    %Done Assuming S =2
    IX = eye(length(A));
    exR = [IX(4*N+1:2:4*N+4*(N-1),:);
           IX(4*N+2:2:4*N+4*(N-1),:)];
    
    res = A + DX\Dy1 + DX\Dy2*kron(eye(2),DV)*kron(eye(2),A)*exR*A;
end

function [DV] = computeDV(X,Para)
    S = 2;
    N = Para.N;
    DV = zeros(2*(N-1),length(X));
    
    beta = Para.beta;
    P = Para.P(1,:);
    U = Para.U;
    i = 1;
    ci = reshape(X(1:(N-1)*S),S,N-1);
    X(1:(N-1)*S) = []; i = i +(N-1)*S; icN = i;
    cN = repmat(reshape(X(1:S),S,1),1,N-1);
    X(1:S) = [];i=i+S;
    X(1:(N-1)*S) = [];i = i +(N-1)*S;
    X(1:S) = [];i=i+S;
    X(1:S*(N-1)) = [];i = i +(N-1)*S;
    X(1:S*(N-1)) = [];i = i +(N-1)*S; imu = i;
    mu = reshape(X(1:(N-1)*S),S,N-1);
    X(1:S*(N-1)) = [];i = i +(N-1)*S; ilambda = i;
    lambda = repmat(reshape(X(1:N-1),1,N-1),S,1);

    [~,uci,~,ucci] = U(ci,0.5*ones(2,N-1));
    [~,ucN,~,uccN] = U(cN,0.5*ones(2,N-1));
    
    Euci = P * uci;
    EucN = repmat(P * ucN,S,1);
    EucN_mu = repmat(P * (mu.*ucN),S,1);
    PP = repmat(P',1,N-1);
    
    
    DVx_dcN = PP.*mu.*uccN./(beta*EucN) - PP.*uccN.*EucN_mu./(beta*EucN.^2);
    DVx_dmu = P.*ucN(:,1)'./(beta*EucN(:,1)'); %1xS
    DVx_dmu = kron(eye(N-1),DVx_dmu); %mu_i only effects Vx_i now have N-1x(N-1)*S
    
    DV(1:N-1,icN:icN+S-1) = DVx_dcN';
    DV(1:N-1,imu:imu+(N-1)*S-1) = DVx_dmu;
    
    DVR_dci = -lambda.*PP.*ucci;
    DVR_dlambda = -Euci; %currently 1xN-1
    DVR_dci = DVR_dci';
    n = 1;
    for i = 1:N-1
        DV(N-1+i, n:i*S) = DVR_dci(i,:);
        n = n + S;
    end
    DV(N:2*N-2,ilambda:ilambda+N-2) = diag(DVR_dlambda);
    

end

function [Dg1 Dg2] = fDerivativetest(f,XSS,YSS)
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

function [DV] = computeDVtest(XSS,Para)
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
        
    DVx_dc2 = mu.*P.*ucc2/(beta*Euc2) - (P.*ucc2*Euc2_mu)/(beta*Euc2^2);
    DVx_dmu = P.*uc2/(beta*Euc2);
    
    DV(1,3:4) = DVx_dc2;
    DV(1,13:14) = DVx_dmu;
    
    
    DVR_dc1 = -lambda*P.*ucc1;
    DVR_dlambda = -Euc1;
    DV(2,1:2) = DVR_dc1;
    DV(2,15) = DVR_dlambda;

end

function [ret] = transform2to3(A2)
    ret = [repmat(A2(1:2,:),2,1);A2(3:4,:);repmat(A2(5:6,:),2,1);A2(7:8,:);repmat(A2(9:10,:),2,1);repmat(A2(11:12,:),2,1);repmat(0.5*A2(13:14,:),2,1);repmat(0.5*A2(15,:),2,1);repmat(0.5*A2(16:17,:),2,1);A2(18:19,:);repmat(A2(20:21,:)*0.5,2,1)];
end