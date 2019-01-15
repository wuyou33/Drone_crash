function calc_lqr(func_LDM,traj_name,P)
% This script calculate lqr gains along the trajectory
traj = load(traj_name);

Q = eyes(length(traj.X{1}));
R = 0.1*eyes(length(traj.U{1}));

A_traj = zeros(length(traj.X{1}),length(traj.X{1}),length(traj.X));
B_traj = zeros(length(traj.X{1}),length(traj.U{1}),length(traj.X));
for i = 1:length(traj.X)
    A_traj(:,:,i) = calc_As(traj.X{i},traj.U{i},P,func_LDM);
    B_traj(:,:,i) = calc_Bs(traj.X{i},traj.U{i},P);
end

end