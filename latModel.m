%% Open-Loop Lateral Dynamic Model for Articulated Tractor Trailer System
function lat = latModel(steer_ang, Vx, dt)

% Author: Tahn Thawainin, AU GAVLAB
%
% Description: Function to simulate open loop dynamics for an
%              articulated tractor trailer system (Ref: Wolfe)
%
% Inputs: steer_ang - wheel steer angle (rad)
%         Vx - longitudinal velocity (m/s)
%         dt - sampling rate
%
% Outputs: lat_model - open loop lateral dynamics data set (SI)

%% Vehicle Parameters
vp = vehParams();

% longitudinal velocity
% Vx = 25;

%% Model Matrix

% revision 4- assumes cos(hitch) = 1 and sin(hitch) = 0
lat.M = [vp.m_t1 + vp.m_t2, -vp.m_t2*(vp.c + vp.d), Vx*(vp.m_t1 + vp.m_t2), -vp.m_t2*vp.d, 0;

     -vp.m_t2*(vp.c + vp.d), vp.j_zz1 + vp.j_zz2 + vp.m_t2*(vp.c + vp.d)^2, ...
     -vp.m_t2*Vx*(vp.c + vp.d), vp.j_zz2 + vp.m_t2*vp.d^2 + vp.m_t2*vp.c*vp.d, 0;

     0, 0, 1, 0, 0;

     -vp.m_t2*vp.d, vp.j_zz2 + vp.m_t2*vp.d^2 + vp.m_t2*vp.c*vp.d, -vp.m_t2*Vx*vp.d, ...
     vp.j_zz2 + vp.m_t2*vp.d^2, 0;

     0, 0, 0, 0, 1];

%% Stiffness Matrix

% matrix elements
k11 = (1/Vx)*(-C1 - C2 - C3 - cos(hitch)*C4 - cos(hitch)*C5);

k12 = (1/Vx)*(-C1*vp.a + C2*vp.b1 + C3*vp.b2 + cos(hitch)*C4*(vp.c + vp.f1*cos(hitch)) ...
       + cos(hitch)*C5*(vp.c + vp.f2*cos(hitch)));

k14 = (1/Vx)*(cos(hitch)^2*vp.f1*C4 + cos(hitch)^2*vp.f1*C5);

k15 = cos(hitch)*C4 + cos(hitch)*C5;

k21 = (1/Vx)*(-C1*vp.a + C2*vp.b1 + C3*vp.b2 + C4*(vp.f1 + vp.c*cos(hitch)) ...
       + C5*(vp.f2 + vp.c*cos(hitch)));

k22 = (1/Vx)*(-C1*vp.a^2 - C2*vp.b1^2 - C3*vp.b2^2 ...
      - (vp.f1 + vp.c*cos(hitch))*C4*(vp.c + vp.f1*cos(hitch)) ...
      - (vp.f2 + vp.c*cos(hitch))*C5*(vp.c + vp.f2*cos(hitch)));

k24 = (1/Vx)*(-(vp.f1 + vp.c*cos(hitch))*C4*vp.f1*cos(hitch)...
      - (vp.f2 + vp.c*cos(hitch))*C5*vp.f2*cos(hitch));

k25 = -(vp.f1 + vp.c*cos(hitch))*C4 - (vp.f2 + vp.c*cos(hitch))*C5;

k41 = (1/Vx)*(vp.f1*C4 + vp.f2*C5);

k42 = (1/Vx)*(-vp.f1*C4*(vp.c + vp.f1*cos(hitch)) - vp.f2*C5*(vp.c + vp.f2*cos(hitch)));

k44 = (1/Vx)*(-vp.f1^2*C4*cos(hitch) - vp.f2^2*C5*cos(hitch));

k45 = -vp.f1*C4 - vp.f2*C5;

% stiffness matrix
lat.K = [k11, k12, 0, k14, k15;
         k21, k22, 0, k24, k25;
         0, 1, 0, 0, 0;
         k41, k42, 0, k44, k45;
         0, 0, 0, 1, 0];

%% Forcing Matrix

lat.F = [cos(hitch)*C1;
         vp.a*cos(hitch)*C1;
         0;
         0;
         0];

%% Continuous State Space Model

% state transition matrix
lat.Ac = M\K;

% input matrix
lat.Bc = M\F;

% observation matrix
lat.Cc = eye(5);

% measurement input matrix
lat.Dc = 0;

% continuos system
lat.sysc = ss(lat.Ac, lat.Bc, lat.Cc, lat.Dc);

%% Discrete State Space Model

% discrete system
lat.sysd = c2d(lat.sysc, dt, 'zoh');

end
