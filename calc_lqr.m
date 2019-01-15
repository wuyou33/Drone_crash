function [K_lqr,A_traj,B_traj] = calc_lqr(X,U,P,func_LDM)
% This script calculate lqr gains along the trajectory specified by X and U

Q = eye(length(X{1}));
R = 0.1*eye(length(U{1}));

K_lqr = zeros(length(U{1}),length(X{1}),length(X));
A_traj = zeros(length(X{1}),length(X{1}),length(X));
B_traj = zeros(length(X{1}),length(U{1}),length(X));
for i = 1:length(X)
    A_traj(:,:,i) = calc_As_5DOF(X{i},U{i},P,func_LDM);
    B_traj(:,:,i) = calc_Bs_5DOF(X{i},U{i},P,func_LDM);
    K_lqr(:,:,i) = lqr(A_traj(:,:,i),B_traj(:,:,i),Q,R);
end

end