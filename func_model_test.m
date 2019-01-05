function LDM = func_model_test(v,a,w)
% user defined longitudinal model to calculate Lift Drag and Moment

D = 0.0*v*v;
L = 0.0*v*v + 1*w;
M = 0.0;

LDM = [L,D,M];
end