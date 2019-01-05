function [X,U] = traj_init_dronecrash(X0,U0,hf,N, P)
% generate init guess of the falling trajectory

dh = hf/N;
g = P(3);

x = zeros(N+1,1); v = zeros(N+1,1); r = zeros(N+1,1); 
a = zeros(N+1,1); q = zeros(N+1,1); w = zeros(N+1,1);

x(1) = X0(1); v(1) = X0(2); r(1) = X0(3); a(1) = X0(4); q(1) = X0(5); w(1) = X0(6);

X = cell(N+1,1); X{1}=X0;
U = cell(N+1,1); U{1}=U0;

% TODO: change to RK4 to improve Numerical stability when gamma and V are
% small

for i = 2:N+1
   x(i) = x(i-1) + dh*(-cot(r(i-1)));
   v(i) = v(i-1) + dh*g/v(i-1);
   r(i) = r(i-1) + dh*g*cot(r(i-1))/v(i-1)^2;
   a(i) = a(i-1);
   q(i) = q(i-1);
   w(i) = w(i-1);
   
   X{i} = [x(i);v(i);r(i);a(i);q(i);w(i)];
   U{i} = U0;
end
end