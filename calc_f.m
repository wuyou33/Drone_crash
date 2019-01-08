function f = calc_f(X,U,P,func_LDM)

v = X(2); 
r = X(3);
a = X(4);
q = X(5);
w = X(6);
u = U(1);

m  = P(1);
Iy = P(2);
g  = P(3);
s  = P(4);

% calculate Lift Drag and Moment from the model
LDM   = func_LDM(v,a,w,q);
L = LDM(1); D = LDM(2); M = LDM(3);

f = zeros(6,1);
% f(x,u):                               % states:
f(1) = -cot(r);                           % X
f(2) = (D+m*g*sin(r))/(m*v*sin(r));       % V
f(3) = -(L-m*g*cos(r))/(m*v^2*sin(r));    % gamma
f(4) = -q/(v*sin(r))-f(3);                % alpha
f(5) = -M/(Iy*v*sin(r));                  % q
f(6) = s/(v*sin(r))*w - s/(v*sin(r))*u;   % omega^2

end