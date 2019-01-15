%% This is the mains cript conducting trajectory optimization 

close all
clear all

addpath(genpath('aero_model'));
addpath('state_equations');
%%
% [As,Bs,bs,Xs,Us,Ps,fs] = calc_derivative();

%% Prameters & Initialization

zf = 30;               % initial height

X0 = [0;                % x0
      3;                % v0
      -30/57.3;         % gamma0
      -30/57.3];        % alpha0
      
U0 =  -30/57.3;         % u0 : alpha_ref

P.m = 0.3;              % mass
P.g = 9.8124;           % g
P.s = 5;                % inner_loop bandwidth
P.Iy = 0.0012;          % moment of inertia

N = 100;
n = length(X0);
m = 1;
dh = zf/N;
xf = 20;

iter_max = 500;
tol = 1e-3;

converg_flag = 0;
thrust_region_scale = 1;
norm_dZ_last = 0;
delta = 1;
theta = 0.8;
nu = 0;
y_last = -inf;
%% aero model
% Func_LDM = @func_model_test;
% Func_LDM = @func_model_18th_Apr;
Func_LDM = @func_model_single_rotor_trim;
%% Init guess
[X,U] =  traj_init_dronecrash(X0,U0,dh,N,P,Func_LDM);
Z = zeros((n+m)*(N+1) + 1,1);
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

% add another variable at the end of Z to regulate xf (chi)
M_calc = [M_calc, zeros(size(M_calc,1),1)];

%% Optimization
% find index of the variables in array Z
id_xf = N*n+1;
id_x  = (0:N)*n+1; id_v  = (0:N)*n+2; id_r  = (0:N)*n+3;
id_a  = (0:N)*n+4;
id_u  = n*(N+1)+(1:N+1);
id_chi = n*(N+1) + m*(N+1) + 1;

% set constraints x+chi >= xf && x-chi <= xf (|x-xf|<=chi)
M_neq = zeros(2,length(Z));
F_neq = zeros(2,1);
M_neq(1,id_xf) = -1; M_neq(1,id_chi) = -1; F_neq(1) = -xf;
M_neq(2,id_xf) =  1; M_neq(2,id_chi) = -1; F_neq(2) =  xf;

for l = 1:100
    % maximum iteration step of X and U
    dx_itr_tol = 1 * thrust_region_scale;
    dv_itr_tol = 1 * thrust_region_scale;
    dr_itr_tol = 3/57.3 * thrust_region_scale;
    da_itr_tol = 3/57.3 * thrust_region_scale;
    du_itr_tol = 3/57.3 * thrust_region_scale;

    % define constraints
    lb    = -inf(size(Z));
    ub    =  inf(size(Z));

    ub(id_x) = min(Z_last(id_x) + dx_itr_tol, inf);
    lb(id_x) = max(Z_last(id_x) - dx_itr_tol,-inf);
    ub(id_v) = min(Z_last(id_v) + dv_itr_tol, inf);
    lb(id_v) = max(Z_last(id_v) - dv_itr_tol, 0);       % V > 0
    ub(id_r) = min(Z_last(id_r) + dr_itr_tol,-5/57.3);  % gamma < -5/57.3
    lb(id_r) = max(Z_last(id_r) - dr_itr_tol,-inf);
    ub(id_a) = min(Z_last(id_a) + da_itr_tol, pi/2);    % alpha < 90
    lb(id_a) = max(Z_last(id_a) - da_itr_tol,-pi/2);    % alpha > -90
    ub(id_u) = min(Z_last(id_u) + du_itr_tol, pi/2);    % alpha_ref < 90
    lb(id_u) = max(Z_last(id_u) - du_itr_tol,-pi/2);    % alpha_ref > -90
    lb(id_chi) = 0;      % chi > 0
    
    % Objective function, 
    f_obj = zeros(size(Z)); 

    % min:1 // max: -1
%     f_obj(id_xf) = -1;
%     f_obj(id_r) = -1;
%     f_obj(id_v) = -1./(Z(id_v).^2.*sin(Z(id_r)));
%     f_obj(id_r) = -cos(Z(id_r))./(Z(id_v).*sin(Z(id_r)).^2);
    f_obj(id_chi) = 1;

    H_obj = zeros(length(Z),length(Z));

    % call the solver
    % [Z,y] = linprog(f_obj,[],[],M_calc,F_calc,lb,ub);
    [Z,y] = quadprog(H_obj,f_obj,M_neq,F_neq,M_calc,F_calc,lb,ub);
    
    norm_dy = norm(y-y_last);
    if nu == -1 && norm_dy < delta/2 * thrust_region_scale^2
        y_last = y; nu = 0;
        break;
    elseif nu == 1 && norm_dy >= delta/2 * thrust_region_scale^2
        y_last = y; nu = 0;
        thrust_region_scale = theta*thrust_region_scale; 
        break;
    end

    if norm_dy > delta/2 * thrust_region_scale^2
        nu = -1;
    elseif norm_dy < delta/2 * thrust_region_scale^2 && norm_dy>tol
        nu = 1;
    else
        nu = 0; converg_flag = 1;
        break;
    end
    
    y_last = y;
    thrust_region_scale = thrust_region_scale * theta^nu;
end

for i = 1:N+1
   X{i} = Z((i-1)*n+1:i*n);
   U{i} = Z(n*(N+1)+i);
end


plot_results(X,U,dh,N,1);

% stop criteria
display([' tol_x = ',num2str(norm(Z-Z_last))]);
display([' tol_y = ',num2str(norm_dy)]);
display([' thrust_region_scale = ',num2str(thrust_region_scale)]);

Z_last = Z;

if num2str(norm(Z-Z_last)) < tol
   converg_flag = 2;
end

if converg_flag == 1 || converg_flag == 2
    display(['Converg! Type = ',num2str(converg_flag)]);
    break;
end
Z(id_chi)
end

%% Plot results
% plot_results(X,U,dh,N,0);
k

return
%% save traj and calculate LQR gains
trim = load('trim_BB2_uva.mat');
for i = 1:N+1
   u_trim = interp2(trim.v,trim.alpha,trim.u',X{i}(2),X{i}(4),'spline');
   X_5D{i} = [X{i};0];
   U_5D{i} = u_trim;
end

K_lqr = calc_lqr(X_5D,U_5D,P,@func_model_single_rotor);
save Traj X_5D U_5D K_lqr