function [As,Bs,bs,Xs,Us,Ps,fs] = calc_derivative()
% This script calculate the derivative of system equation


syms x v r a u m g s
syms D(v,a) L(v,a) M(v,a)


% TODO: add accurate D L M model
% D(v,a,w) = 0.0*v*v;
% L(v,a,w) = 0.0*v*v + 1*w;
% M(v,a,w) = 0.0*D;

%%
% f(x,u):                               % states:
dx = -cot(r);                           % X
dv = (D+m*g*sin(r))/(m*v*sin(r));       % V
dr = -(L-m*g*cos(r))/(m*v^2*sin(r));    % gamma
da = s*a/(v*sin(r))-s*u/(v*sin(r));         % alpha

fs = [dx;dv;dr;da];
Xs = [x;v;r;a];
Us = u;
Ps = [m; g; s];
%%
As = [diff(dx,x), diff(dx,v), diff(dx,r), diff(dx,a);
     diff(dv,x), diff(dv,v), diff(dv,r), diff(dv,a);
     diff(dr,x), diff(dr,v), diff(dr,r), diff(dr,a);
     diff(da,x), diff(da,v), diff(da,r), diff(da,a)];

Bs = [diff(dx,u);
     diff(dv,u);
     diff(dr,u);
     diff(da,u)];

bs = fs - As*Xs - Bs*Us; 
end
