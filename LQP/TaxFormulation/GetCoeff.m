function res = GetCoeff(x,Q,B)
omega_1=x(1);
omega_2=x(2);

Qtautau=Q(1,1);
Qtauc1=Q(1,3);
Qtauc2=Q(1,4);
QTT=Q(2,2);
Qc1c1=Q(3,3);
Qc2c2=Q(4,4);
Btheta_1tau=B(1,2);
Btheta_2tau=B(1,3);
Btheta_1c1=B(3,2);
Btheta_2c2=B(4,3);

eta_1=sqrt(omega_1/Qc1c1)*Qtauc1;
eta_2=sqrt(omega_2/Qc2c2)*Qtauc2;

kappa_1=sqrt(omega_1/Qc1c1)*Btheta_1c1;
kappa_2=sqrt(omega_2/Qc2c2)*Btheta_2c2;

res(1)=omega_1*eta_1*kappa_1-Btheta_1tau;
res(2)=omega_2*eta_2*kappa_2-Btheta_2tau;

end

