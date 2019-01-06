function [Fx, Fz, My] = FM_BB2_longi_18th_Apr(u,w,u1,u2)
% Longitudinal aerodynamic model identified from Bebop2 flight data.  
% u1, u2 are sum of square of front rotor speed and aft rotor speed (rad/s)
% Forces are measured from optitrack.
% This model is identified from Continues Flight Data Set, which is for
% validating the trim curve. 
% Sihao Sun 
% 18th Apr 2018

Fx = u*-2.166864e-01 + u.^2*1.840640e-02 + u.^3* -9.614389e-04 + w*6.170460e-02;
% M = M0 + M1*wsf + M2*wsb;
M0 = u*1.034570e-02 +...
     u.^2* -6.772798e-04+...
     w*8.637144e-03+...
     u.^2.*w*7.173568e-05+...
     u.*w.^2.*2.625090e-04;
M1 = ones(size(u))* 1.520980e-01 +  ...
     u.^2.*1.037535e-03+...
     u.*w *  1.656398e-03 +...
     u * -1.855894e-03;
M2 = ones(size(u))*  -1.630741e-01 + ...
    u * 8.035134e-03 +...
    u.*w*-2.113690e-04 + ...
    u.^2 *  -6.310193e-04;
% T = T1*wsf + T2*wsb
T0 = sign(w).*w.^2 .* 2.975191e-02 + ...
     w.^3 * -3.770790e-03;
T1 = ones(size(u))*  1.670968e+00 + ...
      u.*-8.584644e-02 +...
      u.^2* 2.199807e-03;
T2 = ones(size(u))* 2.145865e+00 + ...
    u.^2 * 1.969960e-02+... 
    w * 7.283182e-02+...
    u.^3 * -6.840546e-04+... 
    u.^3.*w * -1.966397e-04+...
    u.^2.*w *  4.339120e-03;

Fx = Fx;
Fz = -(T0 + T1.*u1/1e6 + T2.*u2/1e6);
My = M0 + M1.*u1/1e6 + M2.*u2/1e6;
My = 0;
end