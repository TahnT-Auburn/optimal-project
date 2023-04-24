%% Lateral Dynamics Simulation for Articulated Tractor Trailer System
%
% Author: Tahn Thawainin, AU GAVLAB
%
% Description: This script calls funtions to simulate the lateral dynamics
%              for an articulated tractor trailer system. Includes an
%              open-loop model and a closed-loop model using internal estimators

clc
clear variables
close all

%% Load Data

% TruckSim data set
ts_data = load('Run107_multi_dblc0.mat');

%% Simulation Specs

% sampling rate (calculated by subtracting TruckSim's event time)
dt = 1/40;

% simulation time
t_sim = ts_data.T_Event;

%% System Inputs from TruckSim 

% average L1 and R1 steer angles(rad)
steer_ang = deg2rad((ts_data.Steer_L1 + ts_data.Steer_R1)/2);

% longitudinal velocity (m/s)
Vx = ts_data.Vx*(1e3/3600);

% lateral velocity (m/s)
Vy = ts_data.Vy*(1e3/3600);

% yaw (rad)
yaw = deg2rad(ts_data.Yaw);

% yaw rate (rad/s)
yaw_rate = deg2rad(ts_data.AVz);

% hitch (rad)
hitch = deg2rad(ts_data.Art_H);

% hitch rated (rad/s)
hitch_rate = deg2rad(ts_data.ArtR_H);

%% Measurements from TruckSim

% measurment noise condition
% 0 - no added noise
% 1 - add measurement noise
meas_noise = 0;

if meas_noise == 0

    % lateral acceleration (m/s^2) - for cornering stiffness RLS
    Ay = 9.81*ts_data.Ay;

elseif meas_noise == 1
    
    % lat accel measurement noise STD
    sigma_Ay = 0.01;

    % lateral acceleration (m/s^2) - for cornering stiffness RLS
    Ay = 9.81*ts_data.Ay + sigma_Ay*randn(1,length(t_sim));

end

%% Tire Model

% TruckSim tire model data set (.csv file)
ts_tiremodel = 'TireFy101.csv';

% vertical load
Fz1 = 5.35493e4/2;
Fz2 = 2.291725e4/2;
Fz3 = 14831/2;
vert_load = [Fz1, Fz2, Fz3];

% call function
ltm = latTireModel(ts_tiremodel, vert_load);

%% Pre-Process Tire Stiffness Estimattion
% Uses clean signals from TruckSim
    
% initialize
x_init = [1e4;
          1e4;
          1e4;
          1e4;
          1e4];

P_init = [10000, 0, 0, 0, 0;...
          0, 10000, 0, 0, 0;...
          0, 0, 10000, 0, 0;...
          0, 0, 0, 1000, 0;...
          0, 0, 0, 0, 1000];

x = x_init;
P = P_init;

% call estimator
for i = 1:length(t_sim)

rls_cs(i) = cornStiff(steer_ang(i), Vx(i), Vy(i), yaw_rate(i), Ay(i),...
                      hitch(i), hitch_rate(i), x, P);

% update states
x = [rls_cs(i).C1; rls_cs(i).C2; rls_cs(i).C3; rls_cs(i).C4; rls_cs(i).C5];

% update covariance
P = [rls_cs(i).P];

end

%% Open-Loop Dynamic Model

% initialize states
Vy = zeros(length(t_sim));
yaw_rate = zeros(length(t_sim));
yaw = zeros(length(t_sim));
hitch_rate = zeros(length(t_sim));
hitch = zeros(length(t_sim));

% initialize RLS CS estimator
x_init = [1e5;
          1e5;
          1e5;
          1e5;
          1e5];

% initial covariance matrix
P_init = [1000, 0, 0, 0, 0;...
          0, 1000, 0, 0, 0;...
          0, 0, 1000, 0, 0;...
          0, 0, 0, 1000, 0;...
          0, 0, 0, 0, 1000];

x = x_init;
P = P_init;

for i = 1:length(t_sim)


end

%% Interface

% display tire model
disp_tire_model = 'true';
% display tire pre-process tire siffness estimator
disp_pre_ts = 'true';


% tire model
if strcmp(disp_tire_model, 'true') == 1
    
    % overlap figure
    % 0 - DONT overlap
    % 1 - overlap
    fig_overlap = 1;

    % TruckSim Tire Model
    figure
    hold on
    plot(ltm.slip_ang, ltm.Fy1, DisplayName='Fz = 7.36e3')
    plot(ltm.slip_ang, ltm.Fy2, DisplayName='Fz = 1.472e4')
    plot(ltm.slip_ang, ltm.Fy3, DisplayName='Fz = 2.943e4')
    plot(ltm.slip_ang, ltm.Fy4, DisplayName='Fz = 4.415e4')
    plot(ltm.slip_ang, ltm.Fy5, DisplayName='Fz = 5.886e4')

    if fig_overlap == 1
    plot(ltm.slip_ang, ltm.Fy_T1,'--', LineWidth = 1.5, DisplayName='Fz_T1')
    plot(ltm.slip_ang, ltm.Fy_T2,'--', LineWidth = 1.5, DisplayName='Fz_T2')
    plot(ltm.slip_ang, ltm.Fy_T3,'--', LineWidth = 1.5, DisplayName='Fz_T3')
    elseif fig_overlap == 0
    end

    hold off
    title('Tire Model')
    xlabel('Slip angle (rad)')
    ylabel('Lateral Force, Fy, (N)')
    zlabel('Verical Load (N)')
    xlim([0.019,pi/2])
    ylim([0, 4.5e4])
%   legend(Location='best')
    grid
    set(gcf, 'color','w')

    % Interpolated tire model
    figure
    hold on
    plot(ltm.slip_ang, ltm.Fy_T1, DisplayName='Fz_T1')
    plot(ltm.slip_ang, ltm.Fy_T2, DisplayName='Fz_T2')
    plot(ltm.slip_ang, ltm.Fy_T3, DisplayName='Fz_T3')
    hold off
    title('System Tire Model')
    xlabel('Slip angle (rad)')
    ylabel('Lateral Force, Fy, (N)')
    zlabel('Verical Load (N)')
    xlim([0.019,pi/2])
    ylim([0, 4.5e4])
%   legend(Location='best')
    grid
    set(gcf, 'color','w')

elseif strcmp(disp_tire_model, 'false') == 1
end

% pre-process tire stiffness estimator
if strcmp(disp_pre_ts, 'true') == 1
    
    figure
    hold on
    plot(t_sim, rls_cs.C1, DisplayName='C1')
    plot(t_sim, rls_cs.C2, DisplayName='C2')
    plot(t_sim, rls_cs.C3, DisplayName='C3')
    plot(t_sim, rls_cs.C4, DisplayName='C4')
    plot(t_sim, rls_cs.C5, DisplayName='C5')
    title('RLS Cornering Stiffness')
    xlabel('Time [s]')
    ylabel('N/rad')
    xlim([2.5,20])
    legend(Location='best')
    grid

elseif strcmp(disp_pre_ts, 'false') == 1
end


