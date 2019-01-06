function check_FM_func_hedging(X,U,dh,LDM_func)
% This function plot the L, D, M of the designated model LDM_func. Also,
% the function check if the model hedging occurs for those identified from
% the flight data.
% Sihao Sun 5th Jan 2019
    
    N = length(X)-1;
    
    t(1) = 0;
    z = 0:dh:dh*N;
    for i = 1:N+1
        [LDM(i,:),flag_hedge(i)] = LDM_func(X{i}(2),X{i}(4),X{i}(6),X{i}(5));
        if i>1
            % estimate time using the trapezoidal rule
            t(i) = t(i-1) - dh/2*(1/(X{i}(2)*sin(X{i}(3)))+1/(X{i-1}(2)*sin(X{i}(3))));
        end
    end
    
    figure(66)
    subplot(2,2,1); plot(z,LDM(:,1)); hold on; ylabel('L [N]'); xlabel('z [m]');
    subplot(2,2,2); plot(z,LDM(:,2)); hold on; ylabel('D [N]'); xlabel('z [m]');
    subplot(2,2,3); plot(z,LDM(:,3)); hold on; ylabel('M [N]'); xlabel('z [m]');
    subplot(2,2,4); plot(z,flag_hedge); hold on; ylabel('hedging'); xlabel('z [m]');

end