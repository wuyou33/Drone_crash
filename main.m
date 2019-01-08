%% This is the mains cript conducting trajectory optimization 

close all
clear all

addpath('aero_model');
% [As,Bs,bs,Xs,Us,Ps,fs] = calc_derivative();

%% Prameters & Initialization

zf = 30;               % initial height

X0 = [0;                % x0
      3;               % v0
      -30/57.3;         % gamma0
      -30/57.3;         % alpha0
      0;                % q0
      0];             % omega0^2 (normalized)

U0 = X0(6);             % u0 (normalized)

P = [0.3;               % mass
     0.0013;            % Iy
     9.8124;            % g
     50];               % act_dyn

umax = 1.2;
umin = 0;

N = 100;
n = length(X0);
m = 1;
dh = zf/N;

iter_max = 100;
tol = 1e-6;
%% aero model
% Func_LDM = @func_model_test;
% Func_LDM = @func_model_18th_Apr;
Func_LDM = @func_model_single_rotor;

%% Init guess
[X,U] =  traj_init_dronecrash(X0,U0,dh,N,P,Func_LDM);
Z = zeros((n+m)*(N+1),1);
for i = 1:N+1
    Z(n*(i-1)+1:n*i) = X{i};
    Z(n*(N+1)+i) = U{i};
end
Z_last = Z;

plot_results(X,U,dh,N,0);
%% Main loop
M = cell(N+1,2*(N+1));
F = cell(N+1,1);
for i = 1:N+1
    for j = 1:N+1
        M{i,j} = zeros(n);
        M{i,N+1+j} = zeros(n,1);
    end
end
M{1,1} = eye(n);
F{1,1} = -2/dh*X0;

for k = 1:iter_max
%% Set matrice for descritization
% M = zeros(n*n,(n+1)*n+(n+1)*m);[]
for i = 2:N+1
    A_i_1   = calc_As(X{i-1},U{i-1},P,Func_LDM);
    A_i     = calc_As(X{i},U{i},P,Func_LDM);
    B_i_1   = calc_Bs(X{i-1},U{i-1},P);
    B_i     = calc_Bs(X{i},U{i},P);
    
    H_i_1   = eye(n) + dh/2*A_i_1;
    H_i     = eye(n) - dh/2*A_i;
    G_i_1   = dh/2*B_i_1;
    G_i     = dh/2*B_i;
    b_i_1   = calc_f(X{i-1},U{i-1},P,Func_LDM) - A_i_1*X{i-1} - B_i_1*U{i-1};
    b_i     = calc_f(X{i},U{i},P,Func_LDM) - A_i*X{i} - B_i*U{i};
    
    M{i,i-1}    = H_i_1;
    M{i,i}      = -H_i;
    M{i,N+i}    = G_i_1;
    M{i,N+1+i}  = G_i;
    
    F{i,1}  =  b_i_1+b_i; 
    
%     display(['i = ', num2str(i)]);
end

M_calc = cell2mat(M);
F_calc = -dh/2*cell2mat(F);

%% Optimization
% find index of the variables in array Z
id_xf = N*n+1;
id_x  = (0:N)*n+1; id_v  = (0:N)*n+2; id_r  = (0:N)*n+3;
id_a  = (0:N)*n+4; id_q  = (0:N)*n+5; id_w  = (0:N)*n+6;
id_u  = n*(N+1)+(1:N+1);

% maximum iteration step of X and U
dx_itr_tol = 1;
dv_itr_tol = 1;
dr_itr_tol = 1/57.3;
da_itr_tol = 1/57.3;
dq_itr_tol = 1/57.3;
dw_itr_tol = 0.01;
du_itr_tol = 0.01;

% define constraints
lb    = -inf((N+1)*(n+m),1);
ub    =  inf((N+1)*(n+m),1);

ub(id_x) = min(Z_last(id_x) + dx_itr_tol, inf);
lb(id_x) = max(Z_last(id_x) - dx_itr_tol,-inf);
ub(id_v) = min(Z_last(id_v) + dv_itr_tol, inf);
lb(id_v) = max(Z_last(id_v) - dv_itr_tol, 0);       % V > 0
ub(id_r) = min(Z_last(id_r) + dr_itr_tol,-5/57.3);  % gamma < -5/57.3
lb(id_r) = max(Z_last(id_r) - dr_itr_tol,-inf);
ub(id_a) = min(Z_last(id_a) + da_itr_tol, pi/2);    % alpha < 90
lb(id_a) = max(Z_last(id_a) - da_itr_tol,-pi/2);    % alpha > -90
ub(id_q) = min(Z_last(id_q) + dq_itr_tol, inf);
lb(id_q) = max(Z_last(id_q) - dq_itr_tol,-inf);
ub(id_w) = min(Z_last(id_w) + dw_itr_tol, umax);    % omega < umax
lb(id_w) = max(Z_last(id_w) - dw_itr_tol, umin);    % omega > umin
ub(id_u) = min(Z_last(id_u) + du_itr_tol, umax);    % u < umax
lb(id_u) = max(Z_last(id_u) - du_itr_tol, umin);    % u > umin

% Objective function, 
f_obj = zeros((N+1)*(n+m),1); 

% min:1 // max: -1
% f_obj(id_xf) = -1;
f_obj(id_u) = -1;

% call the solver
[Z,y] = linprog(f_obj,[],[],M_calc,F_calc,lb,ub);

for i = 1:N+1
   X{i} = Z((i-1)*n+1:i*n);
   U{i} = Z(n*(N+1)+i);
end
plot_results(X,U,dh,N,mod(k,2));

% stop criteria
display([' tol = ',num2str(norm(Z-Z_last))]);
if norm(Z-Z_last) < tol
   break; 
end
Z_last = Z;
y
end

%% Plot results
% plot_results(X,U,dh,N,0);
k