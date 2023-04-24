%% Longitudinal Closed Loop Simulation
function sim = clSim(contrType, PRS, v, p, grade)
% Author: Tahn Thawainin, AU GAVLAB
%
% Description: A function to simulate closed loop longitudinal dynamics and
%              powertrain dynamics
%              
%
% Inputs: contrType - controller type
%                     0) P
%                     1) PD
%                     2) PI
%                     3) PID
%         config - vehicle configuration
%                  0) 3 axle tractor
%                  1) 5 axle unloaded tractor + trailer 
%                  2) 5 axle loaded tractor + trailer
%         PRS - prereference scaling condition
%               0) DONT use prereference scaling
%               1) USE prereference scaling
%         v - current longitudinal velocity
%         p - current vehicle position
%         grade - current grade profile
%
% Ouputs: sim - simulation data set (SI)

% Vehicle Parameters
vp = vehParams();

% Simulation Procedure
pro = lonProc();

% Controller
contr = lonContr(contrType);

%% Longitudinal Forces

% drag force 
sim.F_drag = pro.D*v^2;

% rolling resistance
sim.F_rr = vp.u_rr*vp.m_veh*pro.g*cos(grade);

% grade force
sim.F_grade = vp.m_veh*pro.g*sin(grade);

%% Feedfroward Torque

% feedfoward grade condition - Set false to capture grade influence on
%                              the system
% 0 - false
% 1 - true
ff_grade = 1;

if ff_grade == 1

% feedforward torque
sim.T_ff = (1/pro.scale_factor)*(sim.F_drag + sim.F_rr + sim.F_grade);

elseif ff_grade == 0

% feedforward torque
sim.T_ff = (1/pro.scale_factor)*(sim.F_drag + sim.F_rr);

end

%% Prereference Scaling

if PRS == 0
    v_des = pro.v_des;
elseif PRS == 1
    v_des = (1/contr.DC)*pro.v_des;
end

%% Simulate Controller

% P controller-------------------------------------------------------------
if contrType == 0

% engine torque into plant
sim.T_eng = contr.Kp*(v_des - v) + sim.T_ff;
% sim.T_eng = 2300;

% apply saturation limit
if sim.T_eng > vp.torque_limit_max
    sim.T_eng = vp.torque_limit_max;
elseif sim.T_eng < vp.torque_limit_min
    sim.T_eng = vp.torque_limit_min;
end

% simulate acceleration
sim.accel = (sim.T_eng*(pro.scale_factor) - sim.F_drag - sim.F_rr - ...
            sim.F_grade - pro.B_eff*v)/pro.M_eff;

% simulate velocity
sim.vel = v + sim.accel*pro.dt;

% simulate position
sim.pos = p + sim.vel*pro.dt;

else
end

%% Powertrain Dynamics

% wheel acceleration
sim.wheel_accel = sim.accel/vp.r_eff;

% wheel speed
sim.wheel_speed = sim.vel/vp.r_eff;

% differential speed
sim.diff_speed = sim.wheel_speed;

% gearbox output speed
sim.gearbox_speed = vp.n_d*sim.diff_speed;

% engine speed
sim.engine_speed = vp.n_t(pro.gear)*sim.gearbox_speed;

end