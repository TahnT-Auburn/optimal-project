%% Tractor/Trailer Vehicle Parameters
function vp = vehParams()

% Author: Tahn Thawainin, AU GAVLAB
%
% Description: A function to store vehicle parameters based on vehicle
%              vp.configuration
%
% Inputs: vp.config - vehicle vp.configuration
%                  0) 3 axle tractor
%                  1) 5 axle unloaded tractor + trailer 
%                  2) 5 axle loaded tractor + trailer
%
% Ouputs: vp - vehicle parameter data set (SI)

%% Vehicle vp.configuration
% set vehicle vp.configuration
% 0) 3 axle tractor
% 1) 5 axle unloaded tractor + trailer 
% 2) 5 axle loaded tractor + trailer
vp.config = 1;

%% Powertrain Specs

% transmission gear ratios
vp.n_t = [11.06, 10.2, 7.062, 4.984, 3.966, 2.831, 2.03, 1.47, 1, 0.74];

% differential gear ratio
vp.n_d = 4.4;

% inertias
vp.j_e = 2.75;
vp.j_t = 0.13;
vp.j_ds = 0.012;
vp.j_diff = 0.028;
vp.j_wheel = 1700;

% damping
vp.b_e = 2.21;
vp.b_t = 1.4;
vp.b_diff = 9.7;

% torque limit
vp.torque_limit_max = 3e3;
vp.torque_limit_min = 0;

%% Aerodynamic Specs

% front area
vp.front_area = 7.8016;

% drag coefficient
vp.cd = 0.79;

%% Wheel Specs

% effective wheel radius
vp.r_eff = 0.510;

% rolling resistance coefficient
vp.u_rr = 0.0041; 

%% Mass/Distance/Inertial Specs

% 3 axle tractor-----------------------------------------------------------
if vp.config == 0
% TODO: generate a math model for products of inertia

% tractor mass
vp.m_s1 = 6493;    % sprung mass
vp.m_us_A1 = 570;  % unsprung mass of 1st axel
vp.m_us_A2 = 785;  % unsprung mass of 2nd axel
vp.m_us_A3 = 785;  % unspring mass of 3rd axel
vp.m_t1 = vp.m_s1 + vp.m_us_A1 + vp.m_us_A2 + vp.m_us_A3; % total tractor mass

% total vehicle mass
vp.m_veh = vp.m_t1;

% tractor distances
vp.a = 1.384;      % dist from tractor CG to front axle
vp.b1 = 3.616;     % dist from tractor CG to 1st rear axle
vp.b2 = 4.886;     % dist from tractor CG to 2nd rear axle
vp.c = 4.251;      % dist from tractor CG to hitch point

% tractor radii of gyration
vp.Rx1 = 1.029;
vp.Ry1 = 1.829;
vp.Rz1 = 1.740;

% tractor inertias
vp.j_xx1 = vp.m_s1*vp.Rx1^2;    % roll inertia
vp.j_yy1 = vp.m_s1*vp.Ry1^2;    % pitch inertia
vp.j_zz1 = vp.m_s1*vp.Rz1^2;    % yaw inertia

vp.j_xy1 = 0;       % products of inertia
vp.j_xz1 = 130;
vp.j_yz1 = 0;

% 5 axle unloaded tractor+trailer------------------------------------------
elseif vp.config == 1

% tractor mass
vp.m_s1 = 6493;    % tractor sprung mass
vp.m_us_A1 = 570;  % unsprung mass of 1st axel
vp.m_us_A2 = 785;  % unsprung mass of 2nd axel
vp.m_us_A3 = 785;  % unspring mass of 3rd axel
vp.m_t1 = vp.m_s1 + vp.m_us_A1 + vp.m_us_A2 + vp.m_us_A3; % total tractor mass

% tractor distances
vp.a = 1.384;      % dist from tractor CG to front axle
vp.b1 = 3.616;     % dist from tractor CG to 1st rear axle
vp.b2 = 4.886;     % dist from tractor CG to 2nd rear axle
vp.c = 4.251;      % dist from tractor CG to hitch point

% tractor radii of gyration
vp.Rx1 = 1.029;
vp.Ry1 = 1.829;
vp.Rz1 = 1.740;

% tractor inertias
vp.j_xx1 = vp.m_s1*vp.Rx1^2;    % roll inertia
vp.j_yy1 = vp.m_s1*vp.Ry1^2;    % pitch inertia
vp.j_zz1 = vp.m_s1*vp.Rz1^2;    % yaw inertia

vp.j_xy1 = 0;       % products of inertia
vp.j_xz1 = 130;
vp.j_yz1 = 0;

% trailer mass
vp.m_s2 = 3196;    % trailer sprung mass
vp.m_us_A4 = 665;  % unsprung mass of 4th axel
vp.m_us_A5 = 665;  % unsprung mass of 5th axel
vp.m_t2 = vp.m_s2 + vp.m_us_A4 + vp.m_us_A5; % total trailer mass

% total vehicle mass
vp.m_veh = vp.m_t1 + vp.m_t2;

% trailer distances
vp.d = 7.303;
vp.f1 = 13.124;
vp.f2 = 14.412;

% trailer radii of gyration
vp.Rx2 = 1.765;
vp.Ry2 = 7.322;
vp.Rz2 = 7.505;

% trailer inertias
vp.j_xx2 = vp.m_s2*vp.Rx2^2;    % roll inertia
vp.j_yy2 = vp.m_s2*vp.Ry2^2;    % pitch inertia
vp.j_zz2 = vp.m_s2*vp.Rz2^2;    % yaw inertia

vp.j_xy2 = 0;       % products of inertia
vp.j_xz2 = 0;
vp.j_yz2 = 0;


% 5 axle loaded tractor+trailer--------------------------------------------
elseif vp.config == 2
% TODO: generate math model for different load CG locations in the trailer.
%       This will change the calculation of the inertia properties

% tractor mass
vp.m_s1 = 6493;    % tractor sprung mass
vp.m_us_A1 = 570;  % unsprung mass of 1st axel
vp.m_us_A2 = 785;  % unsprung mass of 2nd axel
vp.m_us_A3 = 785;  % unspring mass of 3rd axel
vp.m_t1 = vp.m_s1 + vp.m_us_A1 + vp.m_us_A2 + vp.m_us_A3; % total tractor mass

% tractor distances
vp.a = 1.384;      % dist from tractor CG to front axle
vp.b1 = 3.616;     % dist from tractor CG to 1st rear axle
vp.b2 = 4.886;     % dist from tractor CG to 2nd rear axle
vp.c = 4.251;      % dist from tractor CG to hitch point

% tractor radii of gyration
vp.Rx1 = 1.029;
vp.Ry1 = 1.829;
vp.Rz1 = 1.740;

% tractor inertias
vp.j_xx1 = vp.m_s1*vp.Rx1^2;    % roll inertia
vp.j_yy1 = vp.m_s1*vp.Ry1^2;    % pitch inertia
vp.j_zz1 = vp.m_s1*vp.Rz1^2;    % yaw inertia

vp.j_xy1 = 0;       % products of inertia
vp.j_xz1 = 130;
vp.j_yz1 = 0;

% load mass
vp.m_l = 17000;

% trailer mass
vp.m_s2 = 3196;    % trailer sprung mass
vp.m_us_A4 = 665;  % unsprung mass of 4th axel
vp.m_us_A5 = 665;  % unsprung mass of 5th axel
vp.m_t2 = vp.m_l + vp.m_s2 + vp.m_us_A4 + vp.m_us_A5; % total trailer mass

% total vehicle mass
vp.m_veh = vp.m_t1 + vp.m_t2;

% trailer distances
vp.d = 7.303;
vp.f1 = 13.124;
vp.f2 = 14.412;

% trailer radii of gyration
vp.Rx2 = 1.765;
vp.Ry2 = 7.322;
vp.Rz2 = 7.505;

% trailer inertias
vp.j_xx2 = vp.m_s2*vp.Rx2^2;    % roll inertia
vp.j_yy2 = vp.m_s2*vp.Ry2^2;    % pitch inertia
vp.j_zz2 = vp.m_s2*vp.Rz2^2;    % yaw inertia

vp.j_xy2 = 0;       % products of inertia
vp.j_xz2 = 0;
vp.j_yz2 = 0;

end

end