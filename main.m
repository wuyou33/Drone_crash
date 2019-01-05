% This is the mains cript conducting trajectory optimization 

% [As,Bs,bs,Xs,Us,Ps,fs] = calc_derivative();

%% Prameters & Initialization

zf = 100;               % initial height

X0 = [0;                % x0
      30;               % v0
      -30/57.3;         % gamma0
      -30/57.3;         % alpha0
      0;                % q0
      2^2];             % omega0^2 (normalized)

U0 = X0(6);             % u0 (normalized)

P = [0.3;               % mass
     0.0013;            % Iy
     9.8124;            % g
     50];               % act_dyn

umax = 10;
umin = 0;

N = 50;
n = length(X0);
m = 1;

iter_max = 1;

% aero model
Func_LDM = @func_model_test;

%% Init guess
[X,U] =  traj_init_dronecrash(X0,U0,zf,N,P);
Z = zeros((n+m)*(N+1),1);
for i = 1:N+1
    Z(n*(i-1)+1:n*i) = X{i};
    Z(n*(N+1)+i) = U{i};
end
%% Main loop
dh = zf/N;
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
    tic
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
id_v  = (0:N)*n+2;
id_r  = (0:N)*n+3;

% define constraints
lb    = -inf((N+1)*(n+m),1);
ub    =  inf((N+1)*(n+m),1);

lb(n*(N+1)+1:end) = umin;
ub(n*(N+1)+1:end) = umax;
ub(id_r)  = -5/57.3;
lb(id_v)  = 0;

% Objective function, 
f_obj = zeros((N+1)*(n+m),1); 

% min:1 // max: -1
f_obj(id_xf) = 1;

% call tje solver
Z = linprog(f_obj,[],[],M_calc,F_calc,lb,ub);

for i = 1:N+1
   X{i} = Z((i-1)*n+1:i*n);
   U{i} = Z(n*(N+1)+i);
end

end

%% Plot results
plot_results(X,U,dh,N,0);
