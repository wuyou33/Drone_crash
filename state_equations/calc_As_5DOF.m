function A = calc_As_5DOF(X,U,P,func_LDM)

v = X(2); 
r = X(3);
a = X(4);
q = X(5);

u = U(1);

m  = P.m;
Iy = P.Iy;
g  = P.g;

% calculate Lift Drag and Moment from the model
LDM   = func_LDM(v,a,u,q);
L = LDM(1); D = LDM(2); M = LDM(3);

% calculate aerodynamic derivative
dv = 0.001; da = 0.001; dq = 0.001;

dLDMdv = (func_LDM(v+dv,a,u,q)-LDM)/dv;
dLDMda = (func_LDM(v,a+da,u,q)-LDM)/da;
dLDMdq = (func_LDM(v,a,u,q+dq)-LDM)/dq;

Lv = dLDMdv(1); Dv = dLDMdv(2); Mv = dLDMdv(3);
La = dLDMda(1); Da = dLDMda(2); Ma = dLDMda(3);
Lq = dLDMdq(1); Dq = dLDMdq(2); Mq = dLDMdq(3);

A = ... 
[ 0,                                                                        0,                                                                 cot(r)^2 + 1,                  0,                                0;
  0,                        Dv/(m*v*sin(r)) - (g*m*sin(r) + D)/(m*v^2*sin(r)),             (g*cos(r))/(v*sin(r)) - (cos(r)*(g*m*sin(r) + D))/(m*v*sin(r)^2),    Da/(m*v*sin(r)),                  Dq/(m*v*sin(r));
  0,                - (2*(g*m*cos(r) - L))/(m*v^3*sin(r)) - Lv/(m*v^2*sin(r)),                         - g/v^2 - (cos(r)*(g*m*cos(r) - L))/(m*v^2*sin(r)^2), -La/(m*v^2*sin(r)),               -Lq/(m*v^2*sin(r));
  0, q/(v^2*sin(r)) + (2*(g*m*cos(r) - L))/(m*v^3*sin(r)) + Lv/(m*v^2*sin(r)), g/v^2 + (q*cos(r))/(v*sin(r)^2) + (cos(r)*(g*m*cos(r) - L))/(m*v^2*sin(r)^2),  La/(m*v^2*sin(r)), Lq/(m*v^2*sin(r)) - 1/(v*sin(r));
  0,                                     M/(Iy*v^2*sin(r)) - Mv/(Iy*v*sin(r)),                                                   (cos(r)*M)/(Iy*v*sin(r)^2),  -Ma/(Iy*v*sin(r)),                -Mq/(Iy*v*sin(r))];    
end