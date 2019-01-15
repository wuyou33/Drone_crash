function f = calc_f(X,U,P,func_LDM)

v = X(2); 
r = X(3);
a = X(4);
u = U(1);

m  = P.m;
g  = P.g;
s  = P.s;

% calculate Lift Drag and Moment from the model
LDM   = func_LDM(v,a);
L = LDM(1); D = LDM(2);

f = zeros(4,1);
% f(x,u):                                 % states:
f(1) = -cot(r);                           % X
f(2) = (D+m*g*sin(r))/(m*v*sin(r));       % V
f(3) = -(L-m*g*cos(r))/(m*v^2*sin(r));    % gamma
f(4) = s/(v*sin(r))*a-s/(v*sin(r))*u;     % alpha

end