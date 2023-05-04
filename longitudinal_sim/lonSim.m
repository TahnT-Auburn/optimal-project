%% Longitudinal Truck Simulation
%
% Author: Tahn Thawainin, AU GAVLAB
%
% Description: Longitudinal simulation environment. This script calls
%              functions to simulate a longitudinal driving scenario
%              
clc
close all
clear variables

%% Vehicle Parameters
% To edit, refer to function {vehParams}
vp = vehParams();

%% Simulation Procedure
% To edit, refer to function {lonProc}
pro = lonProc();

%% Controller
% To edit, refer to function {lonContr}

% set controller type
% 0) P
% 1) PD
% 2) PI
% 3) PID
contrType = 0;
contr = lonContr(contrType);

%% Run Closed Loop Simulation
% To edit, refer to function {clSim}

% initialize
p = pro.pos_init;
v = pro.v_init;

% prereference scaling condition
% 0) DONT prereference scale
% 1) prereference scale
PRS = 1;

for i = 1:length(pro.t_sim)

    % run simulation
    sim(i) = clSim(contrType, PRS, v, p, pro.grade(i));
    
    % update velocity
    v = sim(i).vel;
    
    % update position
    p = sim(i).pos;
    
end

%% Extract Simulation Variables

% position
pos = extractfield(sim,'pos');

% velocity
vel = extractfield(sim,'vel');

% acceleration
accel = extractfield(sim,'accel');

% engine torque
T_eng = extractfield(sim, 'T_eng');

% drag force
F_drag = extractfield(sim, 'F_drag');

% rolling resistance
F_rr = extractfield(sim, 'F_rr');

% grade force
F_grade = extractfield(sim, 'F_grade');

% wheel acceleration
wheel_accel = extractfield(sim, 'wheel_accel');

% wheel speed
wheel_speed = extractfield(sim, 'wheel_speed');

% differential speed
diff_speed = extractfield(sim, 'diff_speed');

% gearbox speed
gearbox_speed = extractfield(sim, 'gearbox_speed');

% engine speed
engine_speed = extractfield(sim, 'engine_speed');
%% Save Data Sets

% save data
save_data = 'true';

if strcmp(save_data, 'true') == 1
    
    % file name
    filename = 'optimal_test1.mat';

    % save file
    save(filename, 'vp','pro','sim')

elseif strcmp(save_data, 'false') == 1
end

%% Interface

% display vehicle info
veh_info = 'true';
% display procedure info
pro_info = 'true';
% display controller info
contr_info = 'true';
% display simulation info
sim_info = 'true';
% display powertrain info
pwrtrain_info = 'true';

% vehicle info
if strcmp(veh_info, 'true') == 1
    
    disp('Vehicle configuration:')

    if vp.config == 0

        disp('3 axle tractor')

        disp('vehicle mass')
        disp(vp.m_veh)

    elseif vp.config == 1

        disp('5 axle unloaded tractor + trailer')
    
        disp('vehicle mass')
        disp(vp.m_veh)

    elseif vp.config == 2

        disp('5 axle loaded tractor + trailer')

        disp('load mass')
        disp(vp.m_l)

        disp('vehicle mass')
        disp(vp.m_veh)
    end

elseif strcmp(veh_info, 'false') == 1
end

% procedure info
if strcmp(pro_info, 'true') == 1

    disp('Procedure Specs:')
    
    disp('sampling rate')
    disp(pro.dt)

    disp('simulation time')
    disp(pro.t_run)
    
    disp('initial position')
    disp(pro.pos_init)

    disp('initial velocity')
    disp(pro.v_init)
    
    disp('desired velocity')
    disp(pro.v_des)
    
    if pro.grade_status == 0
    
    disp('grade status: constant')
   
    disp('grade (deg)')
    disp(pro.grade_deg)

    elseif pro.grade_status == 1
    
    disp('grade status: sinusoidal')

    disp('grade mag (deg)')
    disp(pro.grade_deg)
    
    disp('grade change frequecy')
    disp(pro.grade_freq)
    end


elseif strcmp(pro_info, 'false') == 1
end

% controller info
if strcmp(contr_info, 'true') == 1
    
    disp('Controller Specs:')

    if contrType == 0
    
        disp('Controller Type - P')
        
        disp('Time Constant')
        disp(contr.t_P)

        disp('Controller Gain(s)')
        disp('Kp:')
        disp(contr.Kp)

        disp('DC Gain')
        disp(contr.DC)

    else
    end
    
elseif strcmp(contr_info, 'false') == 1
end

% simulation information
if strcmp(sim_info, 'true') == 1
    
    % sim specs
    figure
    set(gcf,'color','w')
    subplot(4,1,1)
    plot(pro.t_sim, pos)
    title('Position')
    ylabel('m')
    grid

    subplot(4,1,2)
    plot(pro.t_sim, vel)
    title('Velocity')
    ylabel('m/s')
    grid

    subplot(4,1,3)
    plot(pro.t_sim,accel)
    title('Acceleration')
    ylabel('m/s^2')
    grid

    subplot(4,1,4)
    plot(pro.t_sim, T_eng)
    title('Engine Torque')
    ylabel('N-m')
    xlabel('Time [s]')
    grid
    
    % forces
    figure
    set(gcf,'color','w')
    subplot(3,1,1)
    plot(pro.t_sim, F_drag)
    title('Drag Force')
    ylabel('N')
    grid

    subplot(3,1,2)
    plot(pro.t_sim, F_rr)
    title('Rolling Force')
    ylabel('N')
    grid

    subplot(3,1,3)
    plot(pro.t_sim, F_grade)
    title('Grade Force')
    ylabel('N')
    xlabel('Time [s]')
    grid

elseif strcmp(sim_info, 'false') == 1
end

% powertrain information
if strcmp(pwrtrain_info, 'true') == 1
    
    figure
    subplot(3,1,1)
    plot(pro.t_sim, wheel_accel)
    title('Wheel Acceleration')
    ylabel('rad/s^2')
    grid

    subplot(3,1,2)
    plot(pro.t_sim, (30/pi)*wheel_speed)
    title('Wheel Speed')
    ylabel('RPM')
    grid

    subplot(3,1,3)
    plot(pro.t_sim, (30/pi)*engine_speed)
    title('Engine Speed')
    ylabel('RPM')
    xlabel('Time [s]')
    grid
    set(gcf,'color','w')

elseif strcmp(pwrtrain_info, 'false') == 1
end























