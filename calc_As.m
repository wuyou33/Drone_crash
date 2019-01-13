function A = calc_As(X,U,P,func_LDM)

v = X(2); 
r = X(3);
a = X(4);
% q = X(5);
% w = X(6);
% u = U(1);

m  = P(1);
g  = P(2);
s  = P(3);

% calculate Lift Drag and Moment from the model
LDM   = func_LDM(v,a,w,0);
L = LDM(1); D = LDM(2);

% calculate aerodynamic derivative
dv = 0.001; da = 0.001;

dLDMdv = (func_LDM(v+dv,a,w,q)-LDM)/dv;
dLDMda = (func_LDM(v,a+da,w,q)-LDM)/da;

Lv = dLDMdv(1); Dv = dLDMdv(2);
La = dLDMda(1); Da = dLDMda(2);


A = ... 
[ 0,                                                                        0,                                                                 cot(r)^2 + 1,                  0;
  0,                        Dv/(m*v*sin(r)) - (g*m*sin(r) + D)/(m*v^2*sin(r)),             (g*cos(r))/(v*sin(r)) - (cos(r)*(g*m*sin(r) + D))/(m*v*sin(r)^2),    Da/(m*v*sin(r));
  0,                - (2*(g*m*cos(r) - L))/(m*v^3*sin(r)) - Lv/(m*v^2*sin(r)),                         - g/v^2 - (cos(r)*(g*m*cos(r) - L))/(m*v^2*sin(r)^2), -La/(m*v^2*sin(r));
  0,                                  (s*u)/(v^2*sin(r)) - (a*s)/(v^2*sin(r)),                        (s*u*cos(r))/(v*sin(r)^2) - (a*s*cos(r))/(v*sin(r)^2),      s/(v*sin(r))];
 
    
end