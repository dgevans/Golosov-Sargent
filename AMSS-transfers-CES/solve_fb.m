function [n_fb,c_fb,x_fb] =solve_fb(Para,s_)
g=Para.g;
sSize=Para.sSize;
beta=Para.beta;
pi=Para.pi;
der_u_c=Para.der_u_c;
der_u_n=Para.der_u_n;

res_fb=@(n_fb)  der_u_c(n_fb)-der_u_n(n_fb);

n_fb=fsolve(@(n_fb) res_fb(n_fb) , .5*ones(1,sSize));
c_fb=n_fb-g;
uc_fb=der_u_c(n_fb);
un_fb=der_u_n(n_fb);

Euc_fb=sum(pi(s_,:).*uc_fb);
int_fb=uc_fb./(beta*Euc_fb);
x_fb=(uc_fb.*c_fb-un_fb.*n_fb).*(-1+int_fb).^(-1);
end
