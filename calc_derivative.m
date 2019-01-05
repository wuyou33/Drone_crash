function [As,Bs,bs,Xs,Us,Ps,fs] = calc_derivative()
% This script calculate the derivative of system equation


syms x v r a q w u m g s Iy
syms D(v,a,w) L(v,a,w) M(v,a,w)


% TODO: add accurate D L M model
D(v,a,w) = 0.0*v*v;
L(v,a,w) = 0.0*v*v + 1*w;
M(v,a,w) = 0.0*D;

%%
% f(x,u):                               % states:
dx = -cot(r);                           % X
dv = (D+m*g*sin(r))/(m*v*sin(r));       % V
dr = -(L-m*g*cos(r))/(m*v^2*sin(r));    % gamma
da = -q/(v*sin(r))-dr;                  % alpha
dq = -M/(Iy*v*sin(r));                  % q
dw = s/(v*sin(r))*w - s/(v*sin(r))*u;   % omega^2

fs = [dx;dv;dr;da;dq;dw];
Xs = [x;v;r;a;q;w];
Us = u;
Ps = [m; Iy; g; s];
%%
As = [diff(dx,x), diff(dx,v), diff(dx,r), diff(dx,a),diff(dx,q),diff(dx,w);
     diff(dv,x), diff(dv,v), diff(dv,r), diff(dv,a),diff(dv,q),diff(dv,w);
     diff(dr,x), diff(dr,v), diff(dr,r), diff(dr,a),diff(dr,q),diff(dr,w);
     diff(da,x), diff(da,v), diff(da,r), diff(da,a),diff(da,q),diff(da,w);
     diff(dq,x), diff(dq,v), diff(dq,r), diff(dq,a),diff(dq,q),diff(dq,w);
     diff(dw,x), diff(dw,v), diff(dw,r), diff(dw,a),diff(dw,q),diff(dw,w);];

Bs = [diff(dx,u);
     diff(dv,u);
     diff(dr,u);
     diff(da,u);
     diff(dq,u);
     diff(dw,u)];

bs = fs - As*Xs - Bs*Us; 
end
