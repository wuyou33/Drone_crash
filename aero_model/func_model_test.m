function LDM = func_model_test(v,a,w,q)
% user defined longitudinal model to calculate Lift Drag and Moment

D = 0.1*v*v;
L = 0.1*v*v + 0.1*w;
M = 0.01*D - 0.1*L;

LDM = [L,D,M];

end