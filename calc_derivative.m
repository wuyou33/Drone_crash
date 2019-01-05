% This script calculate the derivative of system equation

syms x v r a q w u m g s
syms D(v,a,w) L(v,a,w) M(v,a,w)

%%
% f(x,u):                           % states:
dx = cot(r);                        % X
dv = -D/(m*v*sin(r));               % V
dr = (L-m*g)/(m*v^2*sin(r));        % gamma
da = q/(v*sin(r))-dr;               % alpha
dq = M/(v*sin(r));                  % q
dw = -s*w + s*u;                    % omega^2

%%
A = [diff(dx,x), diff(dx,v), diff(dx,r), diff(dx,a),diff(dx,q),diff(dx,w);
     diff(dv,x), diff(dv,v), diff(dv,r), diff(dv,a),diff(dv,q),diff(dv,w);
     diff(dr,x), diff(dr,v), diff(dr,r), diff(dr,a),diff(dr,q),diff(dr,w);
     diff(da,x), diff(da,v), diff(da,r), diff(da,a),diff(da,q),diff(da,w);
     diff(dq,x), diff(dq,v), diff(dq,r), diff(dq,a),diff(dq,q),diff(dq,w);
     diff(dw,x), diff(dw,v), diff(dw,r), diff(dw,a),diff(dw,q),diff(dw,w);];

B = [0;
     0;
     0;
     0;
     0;
     s];

b = [dx;dv;dr;da;dq;dw] - A*[x;v;r;a;q;w]; 
