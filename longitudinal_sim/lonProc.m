%% Longitudinal Simulation Procedure
function pro = lonProc()
% Author: Tahn Thawainin, AU GAVLAB
%
% Description: A function to specify a longitudinal simulation proedure.
%              Generates simulation specs, drag, grade, gear, and effective terms
%
% Inputs: config - vehicle configuration
%                  0) 3 axle tractor
%                  1) 5 axle unloaded tractor + trailer 
%                  2) 5 axle loaded tractor + trailer              
%
% Ouputs: pro - simulation proedure data set (SI)

% TODO: Add gear shfting and braking procedure

%% Simulation Specs

% sample rate
pro.dt = 0.025;

% run time
pro.t_run = 2;

% simulation time vector
pro.t_sim = 0:pro.dt:pro.t_run;

% initial position
pro.pos_init = 0;

% initial velocity
pro.v_init = 19.44;

% desired velocity
pro.v_des = 19.44;

%% Vehicle Parameters
vp = vehParams();

%% Gravity
pro.g = 9.81;

%% Aerodynamics

% air density
pro.p = 1.206;

% drag term
pro.D = 0.5*pro.p*vp.cd*vp.front_area;

%% Grade

% grade status
% 0) constant grade
% 1) sinusoidal grade
% 2) custom
pro.grade_status = 0;

if pro.grade_status == 0
    
    % input grade (deg)
    pro.grade_deg = 0;

    % road grade (rad)
    pro.grade = deg2rad(pro.grade_deg)*ones(size(pro.t_sim));

elseif pro.grade_status == 1

    % input grade (deg)
    pro.grade_deg = 5;

    % grade change frequency
    pro.grade_freq = 0.25;

    % road grade (rad)
    pro.grade = deg2rad(pro.grade_deg)*sin(pro.grade_freq*(pro.t_sim));

end

%% Gear Shift Schedule

% current gear
pro.gear = 9;

%% Scaling Term

pro.scale_factor = (vp.n_t(pro.gear)*vp.n_d/vp.r_eff);

%% Effective Mass and Damping

% inertial mass
pro.M_i = (vp.j_e*vp.n_t(pro.gear)^2*vp.n_d^2)/vp.r_eff^2 + ...
             (vp.j_t + vp.j_ds)*vp.n_d^2/vp.r_eff^2 + ...
             (vp.j_diff + vp.j_wheel)/vp.r_eff^2;

% effective mass
pro.M_eff = pro.M_i + vp.m_veh;

% effective damping
pro.B_eff = (vp.b_e*vp.n_t(pro.gear)^2*vp.n_d^2)/vp.r_eff^2 + ...
             (vp.b_t*vp.n_d^2/vp.r_eff^2) + ...
             (vp.b_diff/vp.r_eff^2);

end