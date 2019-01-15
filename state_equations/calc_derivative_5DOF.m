function [As,Bs,bs,Xs,Us,Ps,fs] = calc_derivative_5DOF()
% This script calculate the derivative of system equation


syms x v r a q u m g s Iy
syms D(v,a,q,u) L(v,a,q,u) M(v,a,q,u)


% TODO: add accurate D L M model
% D(v,a,u) = 0.0*v*v;
% L(v,a,u) = 0.0*v*v + 1*u;
% M(v,a,u) = 0.0*D;

%%
% f(x,u):                               % states:
dx = -cot(r);                           % X
dv = (D+m*g*sin(r))/(m*v*sin(r));       % V
dr = -(L-m*g*cos(r))/(m*v^2*sin(r));    % gamma
da = -q/(v*sin(r))-dr;                  % alpha
dq = -M/(Iy*v*sin(r));                  % q

fs = [dx;dv;dr;da;dq];
Xs = [x;v;r;a;q];
Us = u;
Ps = [m; Iy; g; s];
%%
As = [diff(dx,x), diff(dx,v), diff(dx,r), diff(dx,a),diff(dx,q);
     diff(dv,x), diff(dv,v), diff(dv,r), diff(dv,a),diff(dv,q);
     diff(dr,x), diff(dr,v), diff(dr,r), diff(dr,a),diff(dr,q);
     diff(da,x), diff(da,v), diff(da,r), diff(da,a),diff(da,q);
     diff(dq,x), diff(dq,v), diff(dq,r), diff(dq,a),diff(dq,q)];

Bs = [diff(dx,u);
     diff(dv,u);
     diff(dr,u);
     diff(da,u);
     diff(dq,u)];

bs = fs - As*Xs - Bs*Us; 
end
