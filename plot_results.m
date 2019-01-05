function plot_results(X,U,dh,N)

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
        t(i) = t(i-1) - dh/2*(1/(v(i)*sin(r(i)))+1/(v(i-1)*sin(r(i-1))));
      end
   end
   
   figure
   subplot(3,2,1);
   plot(z,x); xlabel('Z [m]'); ylabel('x [m]');
   subplot(3,2,2);
   plot(z,v); xlabel('Z [m]'); ylabel('V [m/s]');
   subplot(3,2,3);
   plot(z,r); xlabel('Z [m]'); ylabel('\gamma [rad]');
   subplot(3,2,4);
   plot(z,a); xlabel('Z [m]'); ylabel('\alpha [rad]');
   subplot(3,2,5);
   plot(z,q); xlabel('Z [m]'); ylabel('q [rad/s]');
   subplot(3,2,6);
   plot(z,w); xlabel('Z [m]'); ylabel('omega^2'); hold on; 
   plot(z,u); legend('real','cmd');
   
   figure
   subplot(3,2,1);
   plot(t,x); xlabel('t [s]'); ylabel('x [m]');
   subplot(3,2,2);
   plot(t,v); xlabel('t [s]'); ylabel('V [m/s]');
   subplot(3,2,3);
   plot(t,r); xlabel('t [s]'); ylabel('\gamma [rad]');
   subplot(3,2,4);
   plot(t,a); xlabel('t [s]'); ylabel('\alpha [rad]');
   subplot(3,2,5);
   plot(t,q); xlabel('t [s]'); ylabel('q [rad/s]');
   subplot(3,2,6);
   plot(t,w); xlabel('t [s]'); ylabel('omega^2'); hold on;
   plot(t,u); legend('real','cmd');
   
end