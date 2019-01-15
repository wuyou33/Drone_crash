function [X,U] = traj_init_dronecrash(X0,U0,dh,N, P,func_LDM)
% set U as constant U0 and numerically integrate the initial guess of the
% optimal trajectory.

X = cell(N+1,1); X{1}=X0;
U = cell(N+1,1); 
for i=1:N+1, U{i,1}=U0; end

% TODO: to improve Numerical stability when gamma and V are
% small
for i = 2:N+1
   u = U{i}; 
   X{i} = fsolve(@(x)solve_f(x,u,X{i-1},U{i-1},P,func_LDM,dh),X{i-1}); 
end
end

function f = solve_f(x,u,X,U,P,func_LDM,dh)
f = X + dh/2*(calc_f(X,U,P,func_LDM) + calc_f(x,u,P,func_LDM))-x;
end

function [X_next] =  Rk_4(F,time,X,U,step,param,func_LDM)
    
    e =[0.5,0.5,1.0,1.0,0.5];
	px = X;
	xw = X;
        
	for j=1:4        
		[f] = feval(F,px,U,param,func_LDM);	

        px = xw + e(j)*step*f;
        X = X + e(j+1)*step*f/3.0;        
    end
   
    X_next = X; 
   
end