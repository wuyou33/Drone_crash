function B = calc_Bs_5DOF(X,U,P,func_LDM)

v = X(2); 
r = X(3);
a = X(4);
q = X(5);
u = U(1);
s = P.s;
m = P.m;
Iy= P.Iy;

LDM   = func_LDM(v,a,u,q);

% calculate aerodynamic derivative
du = 0.001;
dLDMdu = (func_LDM(v,a,u+du,q)-LDM)/du;
Lu = dLDMdu(1); Du = dLDMdu(2); Mu = dLDMdu(3);
B =[0;
    Du/(m*v*sin(r));
   -Lu/(m*v*sin(r));
    Lu/(m*v^2*sin(r));
   -Mu/(Iy*v*sin(r))];
end