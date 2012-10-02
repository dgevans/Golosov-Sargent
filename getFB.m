function [c1 c2 l1 l2 y g_y Agent1WageShare]=getFB(Para,s)
psi=Para.psi;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
theta_1=Para.theta_1(s);
theta_2=Para.theta_2(s);
n1=Para.n1;
n2=Para.n2;
g=Para.g;
 options = optimset('Display','off');
 
 l10=.5;
 l20=1-(alpha_2/alpha_1)*(theta_2/theta_1)*l10;
 c10=alpha_1*(n1*theta_1*l10+theta_2*l20-g);
 c20=alpha_2*(n2*theta_1*l10+theta_2*l20-g);
 [x ~, exitflag] =fsolve(@(x) FOCFB(x,s,Para),[c10 c20 l10 l20],options)  ;
 c1=x(1);
 c2=x(2);
 l1=x(3);
 l2=x(4);
 y=l1*n1*theta_1+n2*theta_2*l2;
 g_y=g/y;
 Agent1WageShare=(theta_1*l1)/(theta_2*l2);
end
