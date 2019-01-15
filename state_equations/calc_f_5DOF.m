function f = calc_f_5DOF(X,U,P,func_LDM)

v = X(2); 
r = X(3);
a = X(4);
q = X(5);

u = U(1);

m  = P.m;
Iy = P.Iy;
g  = P.g;
s  = P.s;

% calculate Lift Drag and Moment from the model
LDM   = func_LDM(v,a,u,q);
L = LDM(1); D = LDM(2); M = LDM(3);

f = zeros(5,1);
% f(x,u):                               % states:
f(1) = -cot(r);                           % X
f(2) = (D+m*g*sin(r))/(m*v*sin(r));       % V
f(3) = -(L-m*g*cos(r))/(m*v^2*sin(r));    % gamma
f(4) = -q/(v*sin(r))-f(3);                % alpha
f(5) = -M/(Iy*v*sin(r));                  % q

end