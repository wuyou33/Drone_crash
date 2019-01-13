
function [F,M] = FM_BB2_dml_damage_6DOF(V,pqr,omega,R,l,b,signr)
% This models contains the 3-axis forces + moments of a single rotor and
% the 3-axis forces of a bareframe. The model is not accurate but can be
% used for calculating a preliminary trajectory
% The model has been modified from the archieved model. The pitch moment
% apart from the thrust generated part has been changed to M0. The drag
% model on the airframe is also introduced as M0;
% Sihao Sun 13-Jan-2019

%% load aerodynamic model parameters
persistent model_rotor model_airframe
if isempty(model_rotor)
    model_rotor = load('model_individual_damage_BB2.mat');
    model_airframe = load('Cz_airframe_BB2.mat');
end
%% 
M = model_rotor;
rho = 1.225;

for i=1:4
    if omega(i)<=1
        omega(i) = 0;
    end
end
omega1 = omega(1);
omega2 = omega(2);
omega3 = omega(3);
omega4 = omega(4);

Area = pi*R^2;

u = V(1); v = V(2); w = V(3);

va = sqrt(u^2+v^2+w^2);

if va>=0.01
    alpha = asin(w./va)*57.3;
else
    alpha = 0;
end

%% Calculate aerodynamic normal force of the airframe,

% normal force from the static wind tunnel test
CN_airframe = interp1(model_airframe.AoA_airframe,model_airframe.Cz_airframe,alpha,'spline','extrap');
T0 = -CN_airframe * va * 2;
X0 = -abs(0.1*T0); % assume X0 is proportional to the airframe normal force
M0 = (1.3e-3 + 3e-3*(alpha/57.3) + 1.5e-3*(alpha/57.3)^2)*va^2;

% normal force correctness from the flight data
w_constraint = w;
if w_constraint <=-6
    w_constraint = -6;
elseif w_constraint >2
    w_constraint = 2;
end
u_constraint = u;
v_constraint = v;
uv_scale = 8/sqrt(u^2+v^2);
if uv_scale <=1
    u_constraint = uv_scale*u;
    v_constraint = uv_scale*v;
end
AT0 = P33(u_constraint,abs(v_constraint),w_constraint);
T_corr = AT0*M.k_model{1};

%% Calculate Ct and Cq of each rotor
SL = [1 -1 -1 1];
SM = [1 1 -1 -1];
SN = signr*[1 -1 1 -1];

d = [ l,    -b, 0;   % note that d need to be modified to take rc into account!!
      l,     b, 0;
     -l,     b, 0;
     -l,    -b, 0];
 
V1 = cross(pqr,d(1,:))' + V;
V2 = cross(pqr,d(2,:))' + V;
V3 = cross(pqr,d(3,:))' + V;
V4 = cross(pqr,d(4,:))' + V;

u1 = V1(1); u2 = V2(1); u3 = V3(1); u4 = V4(1);
v1 = V1(2); v2 = V2(2); v3 = V3(2); v4 = V4(2);
w1 = V1(3); w2 = V2(3); w3 = V3(3); w4 = V4(3);

va1 = sqrt(u1^2+v1^2+w1^2);
va2 = sqrt(u2^2+v2^2+w2^2);
va3 = sqrt(u3^2+v3^2+w3^2);
va4 = sqrt(u4^2+v4^2+w4^2); 

vv1 = va1/omega1/R; vv1(isnan(vv1)) = 0; vv1(isinf(abs(vv1))) = 0; vv1(vv1>=0.6) = 0.6;
vv2 = va2/omega2/R; vv2(isnan(vv2)) = 0; vv2(isinf(abs(vv2))) = 0; vv2(vv2>=0.6) = 0.6;
vv3 = va3/omega3/R; vv3(isnan(vv3)) = 0; vv3(isinf(abs(vv3))) = 0; vv3(vv3>=0.6) = 0.6;
vv4 = va4/omega4/R; vv4(isnan(vv4)) = 0; vv4(isinf(abs(vv4))) = 0; vv4(vv4>=0.6) = 0.6;

alpha1 = atan(w1/sqrt(u1^2+v1^2))*57.3; alpha1(isnan(alpha1)) = 0; alpha1(isinf(abs(alpha1))) = 0;
alpha2 = atan(w2/sqrt(u2^2+v2^2))*57.3; alpha2(isnan(alpha2)) = 0; alpha2(isinf(abs(alpha2))) = 0;
alpha3 = atan(w3/sqrt(u3^2+v3^2))*57.3; alpha3(isnan(alpha3)) = 0; alpha3(isinf(abs(alpha3))) = 0;
alpha4 = atan(w4/sqrt(u4^2+v4^2))*57.3; alpha4(isnan(alpha4)) = 0; alpha4(isinf(abs(alpha4))) = 0;

dynhead1 = rho*omega1^2*R^2;
dynhead2 = rho*omega2^2*R^2;
dynhead3 = rho*omega3^2*R^2;
dynhead4 = rho*omega4^2*R^2;

% Ct0 and Cq0 from the static wind tunnel test data
Ct01 = P52(alpha1,vv1)*M.k_Ct0;
Ct02 = P52(alpha2,vv2)*M.k_Ct0;
Ct03 = P52(alpha3,vv3)*M.k_Ct0;
Ct04 = P52(alpha4,vv4)*M.k_Ct0;

Cq01 = P52(alpha1,vv1)*M.k_Cq0;
Cq02 = P52(alpha2,vv2)*M.k_Cq0;
Cq03 = P52(alpha3,vv3)*M.k_Cq0;
Cq04 = P52(alpha4,vv4)*M.k_Cq0;


%% Forces and moments of each rotorM0
T1 = Ct01*dynhead1*Area;
T2 = Ct02*dynhead2*Area;
T3 = Ct03*dynhead3*Area;
T4 = Ct04*dynhead4*Area;

X1 = u1*omega1*M.k_model{3} + SN(1)*v1*omega1*M.k_model{4};
X2 = u2*omega2*M.k_model{3} + SN(2)*v2*omega2*M.k_model{4};
X3 = u3*omega3*M.k_model{3} + SN(3)*v3*omega3*M.k_model{4};
X4 = u4*omega4*M.k_model{3} + SN(4)*v4*omega4*M.k_model{4};

Y1 = SN(1)*u1*omega1*M.k_model{6} + v1*omega1*M.k_model{5};
Y2 = SN(2)*u2*omega2*M.k_model{6} + v2*omega2*M.k_model{5};
Y3 = SN(3)*u3*omega3*M.k_model{6} + v3*omega3*M.k_model{5};
Y4 = SN(4)*u4*omega4*M.k_model{6} + v4*omega4*M.k_model{5};

L1 = SN(1)*u1*omega1*M.k_model{8} + v1*omega1*M.k_model{7};
L2 = SN(2)*u2*omega2*M.k_model{8} + v2*omega2*M.k_model{7};
L3 = SN(3)*u3*omega3*M.k_model{8} + v3*omega3*M.k_model{7};
L4 = SN(4)*u4*omega4*M.k_model{8} + v4*omega4*M.k_model{7};

M1 = u1*omega1*M.k_model{9} + SN(1)*v1*omega1*M.k_model{10};
M2 = u2*omega2*M.k_model{9} + SN(2)*v2*omega2*M.k_model{10};
M3 = u3*omega3*M.k_model{9} + SN(3)*v3*omega3*M.k_model{10};
M4 = u4*omega4*M.k_model{9} + SN(4)*v4*omega4*M.k_model{10};

N1 = SN(1)*(Cq01)*dynhead1*Area;
N2 = SN(2)*(Cq02)*dynhead2*Area;
N3 = SN(3)*(Cq03)*dynhead3*Area;
N4 = SN(4)*(Cq04)*dynhead4*Area;

T = T0+T1+T2+T3+T4;
Fx = X1+X2+X3+X4 + X0;
Fy = Y1+Y2+Y3+Y4;
Mx = SL(1)*b*T1 + SL(2)*b*T2 + SL(3)*b*T3 + SL(4)*b*T4 + L1 + L2 + L3 + L4;
% My = SM(1)*l*T1 + SM(2)*l*T2 + SM(3)*l*T3 + SM(4)*l*T4 + M1 + M2 + M3 + M4;
My = SM(1)*l*T1 + SM(2)*l*T2 + SM(3)*l*T3 + SM(4)*l*T4 + M0;

Mz = b*SL(1)*X1 + b*SL(2)*X2 + b*SL(3)*X3 + b*SL(4)*X4 + ...
     l*SM(1)*Y1 + l*SM(2)*Y2 + l*SM(3)*Y3 + l*SM(4)*Y4 + ...
     N1 + N2 + N3 + N4;
%%
F = [Fx, Fy, -T]';
M = [Mx, My, Mz]';
end
