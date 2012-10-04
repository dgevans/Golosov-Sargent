clc
R=Rold
x=xInit
[x, fvec,~,ifail]=c05qb('BelObjectiveUncondGradNAGBGP',x)
opts = optimset('Algorithm', 'interior-point', 'Display','off','TolX',1e-6);
    xoptguess=x;
    [x, fvec,exitflag]=fminunc(@(x) -Value3cont(x) ,xoptguess,opts)
[x, fvec,~,ifail]=c05qb('BelObjectiveUncondGradNAGBGP',x)

%% GET THE Policy Rules
psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(:,1);
theta_2 = Par.theta(:,2);
g = Par.g;
alpha = Par.alpha;

sigma = 1;
c1_1=x(1);
c1_2=x(2);
c2_1=x(3);

%compute components from unconstrained guess
[c2_2 grad_c2_2] = computeC2_2(c1_1,c1_2,c2_1,R,s_,P,sigma);
[l1 l1grad l2 l2grad] = computeL(c1_1,c1_2,c2_1,c2_2,grad_c2_2,...
    theta_1,theta_2,g,n1,n2);
[btildprime grad_btildprime] = computeBtildeprime(c1_1,c1_2,c2_1,c2_2,grad_c2_2,l1,l2,l1grad,l2grad,...
    u2btild,s_,psi,beta,P);

% x' - definition
u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime


R=Rold
f(x)

f1= @(y) sum((f([y x(2) x(3)]).^2))^.5
ygrid=linspace(x(1)*.9,x(1)*1.1,100)

for yind=1:100
    res(yind)=f1(ygrid(yind));
end

subplot(3,1,1)
plot(ygrid,res)

f2= @(y) sum((f([x(1) y x(3)])).^2)^.5
ygrid=linspace(x(2)*.9,x(2)*1.1,100)

for yind=1:100
    res(yind)=f2(ygrid(yind));
end

subplot(3,1,2)
plot(ygrid,res)


f3= @(y) sum((f([x(1) x(2) y])).^2)^.5
ygrid=linspace(x(3)*.9,x(3)*1.1,100)

for yind=1:100
    res(yind)=f3(ygrid(yind));
end

subplot(3,1,3)
plot(ygrid,res)