function [LDM,flag_hedge] = func_model_18th_Apr(v,a,w,q)
% user defined longitudinal model to calculate Lift Drag and Moment

vx = v*cos(a);
vz = v*sin(a);
% u = sqrt(abs(w))*1000;
u = 2*(w*1000)^2;

% hedge the vx and vz to where the following model supports
flag_hedge = 1;
if vx > 15, vx = 15; elseif vx < -15, vx = -15; else flag_hedge=0; end
if vz > 5,  vz = 5;  elseif vz<-5, vz = -5; else flag_hedge=0; end

[Fx, Fz, My] = FM_BB2_longi_18th_Apr(vx,vz,0,u);

D = -Fx*cos(a) - Fz*sin(a);
L =  Fx*sin(a) - Fz*cos(a); 
M = My;

LDM = [L,D,M];

end