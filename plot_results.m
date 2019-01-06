function plot_results(X,U,dh,N,cls)

   z = 0:dh:dh*N;
   x = zeros(N+1,1);
   v = zeros(N+1,1);
   r = zeros(N+1,1);
   a = zeros(N+1,1);
   q = zeros(N+1,1);
   w = zeros(N+1,1);
   u = zeros(N+1,1);
   t = zeros(N+1,1);
   for i = 1:N+1
      x(i) = X{i}(1); 
      v(i) = X{i}(2);
      r(i) = X{i}(3);
      a(i) = X{i}(4);
      q(i) = X{i}(5);
      w(i) = X{i}(6);
      u(i) = U{i}(1);
      if i>1
        % estimate time using the trapezoidal rule
        t(i) = t(i-1) - dh/2*(1/(v(i)*sin(r(i)))+1/(v(i-1)*sin(r(i-1))));
      end
   end
   
   if cls == 2
       close all;
   end
   
   figure(60)
   subplot(3,2,1);
   plot(z,x); xlabel('Z [m]'); ylabel('x [m]'); hold on; if cls == 1, hold off; end
   subplot(3,2,2);
   plot(z,v); xlabel('Z [m]'); ylabel('V [m/s]'); hold on; if cls == 1, hold off; end
   subplot(3,2,3);
   plot(z,r); xlabel('Z [m]'); ylabel('\gamma [rad]'); hold on; if cls == 1, hold off; end
   subplot(3,2,4);
   plot(z,a); xlabel('Z [m]'); ylabel('\alpha [rad]'); hold on; if cls == 1, hold off; end
   subplot(3,2,5);
   plot(z,q); xlabel('Z [m]'); ylabel('q [rad/s]'); hold on; if cls == 1, hold off; end
   subplot(3,2,6);
   plot(z,w); hold on; if cls == 1, hold off; end
   xlabel('Z [m]'); ylabel('omega^2'); 
   plot(z,u); hold on; 
   legend('real','cmd');  if cls == 1, hold off; end
   
%    figure(61)
%    subplot(3,2,1);hold on; 
%    plot(t,x); xlabel('t [s]'); ylabel('x [m]'); hold on;
%    subplot(3,2,2);
%    plot(t,v); xlabel('t [s]'); ylabel('V [m/s]'); hold on;
%    subplot(3,2,3);
%    plot(t,r); xlabel('t [s]'); ylabel('\gamma [rad]'); hold on;
%    subplot(3,2,4);
%    plot(t,a); xlabel('t [s]'); ylabel('\alpha [rad]'); hold on;
%    subplot(3,2,5);
%    plot(t,q); xlabel('t [s]'); ylabel('q [rad/s]'); hold on;
%    subplot(3,2,6);
%    plot(t,w); xlabel('t [s]'); ylabel('omega^2'); hold on;
%    plot(t,u); legend('real','cmd');
   
end