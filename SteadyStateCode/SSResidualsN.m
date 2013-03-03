function [ res ] = SSResidualsN( X,Para )
%UNTITLED Summary of this function goes here
%   X will be a ((4S+2)N-S-4,1)


%general setup  variables are SxN-1
P = Para.P(1,:);
S = length(P);
N = Para.N;
beta = Para.beta;
alpha = Para.alpha(:)';%ensure row vector
alpha_i = repmat(Para.alpha(1:N-1),S,1);
alpha_N = alpha(N);
if numel(Para.theta) == N
    theta = Para.theta(:)';%ensure row vector
    theta_i = repmat(theta(1:N-1),S,1);
    theta_N = repmat(theta(N),S,N-1);
else
    theta_i = theta(:,1:N-1);
    theta_N = repmat(theta(:,N),1,N-1);
end
n = Para.n(:)';
n_i = repmat(n(1:N-1),S,1);
n_N = n(N);

g = Para.g(:);
U = Para.U;

ci = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];
cN = repmat(reshape(X(1:S),S,1),1,N-1);
X(1:S) = [];
li = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];
lN = repmat(reshape(X(1:S),S,1),1,N-1);
X(1:S) = [];
x = repmat(reshape(X(1:N-1),1,N-1),S,1);
X(1:N-1) = [];
R = repmat(reshape(X(1:N-1),1,N-1),S,1);
X(1:N-1) = [];
mu = repmat(reshape(X(1:N-1),1,N-1),S,1);
X(1:N-1) = [];
lambda = repmat(reshape(X(1:N-1),1,N-1),S,1);
X(1:N-1) = [];
phi = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];
xi = repmat(reshape(X(1:S),S,1),1,N-1);
X(1:S) = [];
rho = reshape(X(1:(N-1)*S),S,N-1);
X(1:(N-1)*S) = [];

[~,uci,uli,ucci,ulli] = U(ci,li);
[~,ucN,ulN,uccN,ullN] = U(cN,lN);

Euci = repmat(P * uci,S,1);
EucN = repmat(P * ucN,S,1); %Constant so I don't need to transform

res = [];

resTemp = x.*ucN./(beta*EucN) - ucN.*(cN-ci) - x + R.*li.*uli - lN.*ulN;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = theta_N.*R.*uli - theta_i.*ulN;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = sum(theta_i.*li-ci,2) + theta_N(:,1).*(lN(:,1)-cN(:,1)) - g;
res = [res;resTemp];

resTemp = R.*uci-ucN;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = alpha_i.*uci + mu.*ucN-xi.*n_i + rho.*R.*ucci;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = alpha_N(:,1)*ucN(:,1) - sum(mu.*(uccN.*(cN-ci)+ucN),2) - n_N*xi(:,1) - sum(rho.*uccN,2);
res = [res;resTemp];

resTemp = alpha_i.*uli + mu.*R.*( uli + li.*ulli ) + phi.*theta_N.*R.*ulli + xi.*theta_i.*n_i;
res = [res;reshape(resTemp,(N-1)*S,1)];

resTemp = alpha_N(:,1).*ulN(:,1) - sum( mu.*(ulN + lN.*ullN) ,2) - sum(phi.*theta_i.*ullN,2)...
    + xi(:,1).*theta_N(:,1).*n_N(:,1);
res = [res;resTemp];

resTemp = -beta.*lambda.*Euci + mu.*li.*uli + lambda.*uci + phi.*uli.*theta_N + rho.*uci;
res = [res;reshape(resTemp,(N-1)*S,1)];



end

