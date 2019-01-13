function [LDM,flag_hedge] = func_model_single_rotor_trim(v,a)
% user defined longitudinal model to calculate Lift and Drag in the trim
% condition.

flag_hedge = 0;

trim = load('trim_BB2_uva.mat');

u = interp2(trim.v,trim.alpha,trim.u',v,a,'extrap')*1000;
vx = v*cos(a);
vz = v*sin(a);

[F,M] = FM_BB2_dml_damage_6DOF([vx,0,vz]',[0,0,0]',[0,0,u,u]',0.075,0.1150,0.0875,-1);
Fx = F(1); Fz = F(3);
My = M(2);

D = -Fx*cos(a) - Fz*sin(a);
L =  Fx*sin(a) - Fz*cos(a); 
M = My;

LDM = [L,D,M];

end