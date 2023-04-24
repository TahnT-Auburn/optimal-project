%% RLS Cornering Stiffness Estimator with Forgetting Factor 
function rls_cs = cornStiff(steer_ang, Vx, Vy, r, Ay, hitch, hitch_rate, x, P)
% Author: 
%           Tahn Thawainin, AU GAVLAB
%
% Description: 
%           A function to estimate tire cornering stiffnesses for a
%           5-axle tractor using a LS/RLS with forgetting factor
%
% Inputs: 
%           steer_ang - steer angle at front axle (rad)
%           VH - longitudinal velocity (m/s)
%           Vy - lateral velocity (m/s)
%           r - yaw rate (rad/s)
%           Ay - lateral acceleration (m/s^2)
%           hitch - hitch angle (rad)
%           hitch_rate - hitch rate (rad/s)
%           phi_init - initalized state matriHconfig
%           P_init - initalized covariance matriH
%
% Outputs: 
%           rls_cs - Cornering stiffness estimator data set (SI)

%% Vehicle Parameters
vp = vehParams();

%% RLS

% forgetting factor
lambda = 0.999;

% measurements
y = Ay;

% scaling vector
H = (1/(vp.m_t1 + vp.m_t2)).*[steer_ang*cos(steer_ang) - cos(steer_ang)*((Vy + vp.a*r)/Vx),...
                              (-Vy + vp.b1*r)/Vx,...
                              (-Vy + vp.b2*r)/Vx,...
                              ((-cos(hitch)*Vy + cos(hitch)*(vp.c + vp.f1*cos(hitch))*r + cos(hitch)^2*vp.f1*hitch_rate)/Vx) + hitch*cos(hitch), ...
                              ((-cos(hitch)*Vy + cos(hitch)*(vp.c + vp.f2*cos(hitch))*r + cos(hitch)^2*vp.f2*hitch_rate)/Vx) + hitch*cos(hitch)];

% update gain matrix
rls_cs.L = P*H'/(lambda + H*P*H');

% update states
rls_cs.innov = (y - H*x);
x = x + rls_cs.L*(rls_cs.innov);

% Update covariance matrix
rls_cs.P = (eye(5) - rls_cs.L*H)*P*(1/lambda);

% %     direct least squares
%     innov = 1;
%     while abs(innov) > 0.1
% 
% %     update gain matrix
%     L = P*H'./(lambda + H*P*H');
% 
% %     update states
%     innov = (y - H*phi);
%     x= x+ pinv(H)*innov;
%     
% %     Update covariance matrix
%     P = (eye(5) - L*H)*P*(1/lambda);

%     end

% Cornering stiffnesses 
rls_cs.C1 = x(1); 
rls_cs.C2 = x(2);
rls_cs.C3 = x(3);
rls_cs.C4 = x(4);
rls_cs.C5 = x(5);

end