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

%% Vehicle Parameters
vp = vehParams();

%% Simulation Specs

% sampling rate (calculated by subtracting TruckSim's event time)
dt = 1/40;

% simulation time
t_sim = ts_data.T_Event;

%% Signals from TruckSim 

% average L1 and R1 steer angles(rad)
steer_ang = deg2rad((ts_data.Steer_L1 + ts_data.Steer_R1)/2);
steer_ang(1:100) = 0;

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

% X global position
Xo = ts_data.Xo;

% Y global position
Yo = ts_data.Yo;

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

% TruckSim axle slip angles (rad)
sa1 = -deg2rad((ts_data.AlphaL1i + ts_data.AlphaR1i)./2);
sa2 = -deg2rad((ts_data.AlphaL2i + ts_data.AlphaR2i)./2);
sa3 = -deg2rad((ts_data.AlphaL3i + ts_data.AlphaR3i)./2);
sa4 = -deg2rad((ts_data.AlphaL4i + ts_data.AlphaR4i)./2);
sa5 = -deg2rad((ts_data.AlphaL5i + ts_data.AlphaR5i)./2);

% truncated slip angles (rad)
sa23 = (sa2 + sa3)./2;
sa45 = (sa4 + sa5)./2;

%% Pre-Process Tire Stiffness Estimattion
% Uses clean signals from TruckSim

% truncate axle condition
% 0 - DONT truncate axles
% 1 - turncate axles
trunc_axle = 0;

% call estimator
rls_cs = cornStiff(ts_data, trunc_axle);

if trunc_axle == 0
% extractfields

C1_est = extractfield(rls_cs, 'C1');    % estimated stiffnesses
C2_est = extractfield(rls_cs, 'C2');
C3_est = extractfield(rls_cs, 'C3');
C4_est = extractfield(rls_cs, 'C4');
C5_est = extractfield(rls_cs, 'C5');

sa1_est = extractfield(rls_cs, 'sa1');  % estimated slip angle
sa2_est = extractfield(rls_cs, 'sa2');
sa3_est = extractfield(rls_cs, 'sa3');
sa4_est = extractfield(rls_cs, 'sa4');
sa5_est = extractfield(rls_cs, 'sa5');

elseif trunc_axle == 1
% extractfields
C1_est = extractfield(rls_cs, 'C1');
C2_est = extractfield(rls_cs, 'C2');
C3_est = extractfield(rls_cs, 'C3');

sa1_est = extractfield(rls_cs, 'sa1');
sa2_est = extractfield(rls_cs, 'sa2');
sa3_est = extractfield(rls_cs, 'sa3');

end
%% Open-Loop Dynamic Model Simulation

% initialize states
Vy_est_ol = 0;
yaw_rate_est_ol = 0;
yaw_est_ol = 0;
hitch_rate_est_ol = 0;
hitch_est_ol = 0;

% initialize position
pos_X_ol = 0;
pos_Y_ol = 0;

x_init = [Vy_est_ol;
          yaw_rate_est_ol;
          yaw_est_ol;
          hitch_rate_est_ol;
          hitch_est_ol];

x = x_init;

% tire stiffness values
CS = [ltm.A1_cs, ltm.A23_cs, ltm.A23_cs, ltm.A45_cs, ltm.A45_cs]; % truth?
% CS = [4.536622535892314e5, 4.826699979105655e5,4.826699979105655e5,6.695921441312158e4,6.695921441312158e4];
Vx_const = 19.4;

for i = 1:length(t_sim)

% input
u = steer_ang(i);

% update dynamics
lat_ol(i) = latModel(steer_ang(i), Vx(i), hitch_est_ol(i), dt, CS);

% simulate dynamics
xd = lat_ol(i).sysc.A*x + lat_ol(i).sysc.B*u;

% update states
x = x + xd.*dt;

% position integration
pos_X_ol = pos_X_ol + (Vx(i)*cos(x(3)) - x(1)*sin(x(3)))*dt;
pos_Y_ol = pos_Y_ol + (Vx(i)*sin(x(3)) + x(1)*cos(x(3)))*dt;

% siphon variables
Vy_est_ol(i+1,:) = x(1);
yaw_rate_est_ol(i+1,:) = x(2);
yaw_est_ol(i+1,:) = x(3);
hitch_rate_est_ol(i+1,:) = x(4);
hitch_est_ol(i+1,:) = x(5);

pos_X_ol_est(i+1,:) = pos_X_ol;
pos_Y_ol_est(i+1, :) = pos_Y_ol;

end

%% Measurements from TruckSim

% measurment noise condition
% 0 - no added noise
% 1 - add measurement noise
meas_noise = 0;

if meas_noise == 0

    % lateral acceleration (m/s^2)
    Ay_meas = 9.81*ts_data.Ay;    
    % yaw rate
    yaw_rate_meas = yaw_rate;

elseif meas_noise == 1
    
    % lat accel measurement noise STD
    sigma_Ay = 0.01;
    % yaw rate measurement noise STD
    sigma_yr = 0.01;

    % lateral acceleration (m/s^2) 
    Ay_meas = 9.81*ts_data.Ay + sigma_Ay*randn(1,length(t_sim));
    % yaw rate
    yaw_rate_meas = yaw_rate + sigma_yr*randn(1,length(t_sim));
end

%% Closed Loop Kalman Filter

% PROCESS AND MEASUREMENT NOISE--------------------------------------------

Q = [0.0001, 0, 0, 0, 0;
     0, 0.1, 0, 0, 0;
     0, 0, 0.1, 0, 0;
     0, 0, 0, 0.1, 0;
     0, 0, 0, 0, 0.001];

R = [0.001, 0;
     0, 0.001];

% INITIALIZE---------------------------------------------------------------

% initialize states
Vy_est_cl = 0;
yaw_rate_est_cl = 0;
yaw_est_cl = 0;
hitch_rate_est_cl = 0;
hitch_est_cl = 0;

% initialize position
pos_X_cl = 0;
pos_Y_cl = 0;

x_init = [Vy_est_cl;
          yaw_rate_est_cl;
          yaw_est_cl;
          hitch_rate_est_cl;
          hitch_est_cl];

P_init = [1, 0, 0, 0, 0;...
          0, 1, 0, 0, 0;...
          0, 0, 1, 0, 0;...
          0, 0, 0, 1, 0;...
          0, 0, 0, 0, 1];

x = x_init;
P = P_init;

% TIME UPDATE--------------------------------------------------------------

for i = 1:length(t_sim)

% input
u = steer_ang(i);

% measurements
y = [Ay_meas(i);
     yaw_rate_meas(i)];

% update dynamics
lat_ol(i) = latModel(steer_ang(i), Vx(i), hitch_est_cl(i), dt, CS);

% simulate dynamics
x = lat_ol(i).sysd.A*x + lat_ol(i).sysd.B*u;

% update states
% x = x + xd.*dt;

% priori covariance
P = lat_ol(i).sysd.A*P*lat_ol(i).sysd.A' + Q;

% MEASUREMENT UPDATE-------------------------------------------------------

% observation matrix
H = [0, Vx(i), 0, 0, 0;
     0, 1, 0, 0, 0];

% kalman gain
K = P*H'/(H*P*H' + R);

% posteriori covariance
P = (eye(5) - K*H)*P;

% measurement prediction innovation
innov = (y - H*x);

% state correction update
x = x + K*(innov);

% SIPHON VARIABLES---------------------------------------------------------

% position integration
pos_X_cl = pos_X_cl + (Vx(i)*cos(x(3)) - x(1)*sin(x(3)))*dt;
pos_Y_cl = pos_Y_cl + (Vx(i)*sin(x(3)) + x(1)*cos(x(3)))*dt;

% siphon variables
Vy_est_cl(i+1,:) = x(1);
yaw_rate_est_cl(i+1,:) = x(2);
yaw_est_cl(i+1,:) = x(3);
hitch_rate_est_cl(i+1,:) = x(4);
hitch_est_cl(i+1,:) = x(5);

pos_X_cl_est(i+1,:) = pos_X_cl;
pos_Y_cl_est(i+1, :) = pos_Y_cl;

end

%% Interface

% display tire model
disp_tire_model = 'false';
% display tire pre-process tire siffness estimator
disp_pre_ts = 'false';
% display true states
disp_true_states = 'false';
% dispaly open loop dynamics
disp_ol = 'true';
% display closed loop dynamics
disp_cl = 'true';
% dispaly position solition
disp_pos = 'true';

% tire model---------------------------------------------------------------
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

% pre-process tire stiffness estimator-------------------------------------
if strcmp(disp_pre_ts, 'true') == 1
    
    % ALL axles
    if trunc_axle == 0
    
    % estimated tire stiffness
    figure
    hold on
    plot(t_sim, C1_est, DisplayName='C1')
    plot(t_sim, C2_est, DisplayName='C2')
    plot(t_sim, C3_est, DisplayName='C3')
    plot(t_sim, C4_est, DisplayName='C4')
    plot(t_sim, C5_est, DisplayName='C4')
    title('RLS Cornering Stiffness')
    xlabel('Time [s]')
    ylabel('N/rad')
    xlim([2.5,20])
    legend(Location='best')
    grid
    set(gcf, 'color','w')
    
    % Lateral Force comparison (Fy1)
    figure
    hold on
    plot(t_sim, ts_data.Fy_A1, DisplayName='TS')
    plot(t_sim, ltm.A1_cs*sa1, DisplayName='Sim Truth')
%     plot(t_sim, C1_est(end)*sa1_est, DisplayName='RLS')
    plot(t_sim, ltm.A1_cs*sa1_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 1')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')
    
    % Lateral Force comparison (Fy2)
    figure
    hold on
    plot(t_sim, ts_data.Fy_A2, DisplayName='TS')
    plot(t_sim, ltm.A23_cs*sa2, DisplayName='Sim Truth')
%     plot(t_sim, C2_est(end)*sa2_est, DisplayName='RLS')
    plot(t_sim, ltm.A23_cs*sa2_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 2')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')

    % Lateral Force comparison (Fy3)
    figure
    hold on
    plot(t_sim, ts_data.Fy_A3, DisplayName='TS')
    plot(t_sim, ltm.A23_cs*sa3, DisplayName='Sim Truth')
%     plot(t_sim, C3_est(end)*sa3_est, DisplayName='RLS')
    plot(t_sim, ltm.A23_cs*sa3_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 3')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')

    % Lateral Force comparison
    figure
    hold on
    plot(t_sim, ts_data.Fy_A4, DisplayName='TS')
    plot(t_sim, ltm.A45_cs*sa4, DisplayName='Sim Truth')
%     plot(t_sim, C4_est(end)*sa4_est, DisplayName='RLS')
    plot(t_sim, ltm.A45_cs*sa4_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 4')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')
    
    figure
    hold on
    plot(t_sim, ts_data.Fy_A5, DisplayName='TS')
    plot(t_sim, ltm.A45_cs*sa5, DisplayName='Sim Truth')
%     plot(t_sim, C5_est(end)*sa5_est, DisplayName='RLS')
    plot(t_sim, ltm.A45_cs*sa5_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 5')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')

    % truncated axles
    elseif trunc_axle == 1
    figure
    hold on
    plot(t_sim, C1_est, DisplayName='C1')
    plot(t_sim, C2_est, DisplayName='C2')
    plot(t_sim, C3_est, DisplayName='C3')
    title('RLS Cornering Stiffness')
    xlabel('Time [s]')
    ylabel('N/rad')
    xlim([2.5,20])
    legend(Location='best')
    grid
    set(gcf, 'color','w')

    end

elseif strcmp(disp_pre_ts, 'false') == 1
end

% true states--------------------------------------------------------------
if strcmp(disp_true_states, 'true') == 1

    figure
    subplot(5,1,1)
    plot(t_sim, Vy, DisplayName='Vy')
    title('True Lateral Velocity')
    ylabel('m/s')
    grid

    subplot(5,1,2)
    plot(t_sim, yaw_rate, DisplayName='yaw rate')
    title('True Yaw Rate')
    ylabel('rad/s')
    grid

    subplot(5,1,3)
    plot(t_sim, yaw, DisplayName='yaw')
    title('True Yaw')
    ylabel('rad')
    grid

    subplot(5,1,4)
    plot(t_sim, hitch_rate, DisplayName='hitch rate')
    title('True Hitch Rate')
    ylabel('rad/s')
    grid

    subplot(5,1,5)
    plot(t_sim, hitch, DisplayName='hitch')
    title('True Hitch')
    xlabel('Time [s]')
    ylabel('rad')
    grid
    set(gcf, 'color','w')

elseif strcmp(disp_true_states, 'false') == 1
end

% open loop----------------------------------------------------------------
if strcmp(disp_ol, 'true') == 1

    figure
    subplot(5,1,1)
    hold on
    plot(t_sim, Vy, DisplayName='TruckSim')
    plot(t_sim, Vy_est_ol(1:end-1), DisplayName='Model')
    hold off
    title('Tractor Lateral Velocity $V_{y}$', 'Interpreter','latex')
    ylabel('m/s')
    legend(Location='best')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)
    
    subplot(5,1,2)
    hold on
    plot(t_sim, rad2deg(yaw_rate), DisplayName='TruckSim')
    plot(t_sim, rad2deg(yaw_rate_est_ol(1:end-1)), DisplayName='Model')
    hold off
    title('Yaw Rate $\dot{\psi}$', 'Interpreter', 'latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,3)
    hold on
    plot(t_sim, rad2deg(yaw), DisplayName='TruckSim')
    plot(t_sim, rad2deg(yaw_est_ol(1:end-1)), DisplayName='Model')
    hold off
    title('Yaw $\psi$', 'Interpreter','latex')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)
   
    subplot(5,1,4)
    hold on
    plot(t_sim, rad2deg(hitch_rate), DisplayName='hitch rate')
    plot(t_sim, rad2deg(hitch_rate_est_ol(1:end-1)), DisplayName='hitch rate')
    hold off
    title('Hitch Rate $\dot{\gamma}$', 'Interpreter','latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,5)
    hold on
    plot(t_sim, rad2deg(hitch), DisplayName='hitch')
    plot(t_sim, rad2deg(hitch_est_ol(1:end-1)), DisplayName='hitch')
    hold off
    title('Hitch $\gamma$', 'Interpreter','latex')
    xlabel('Time [s]')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

elseif strcmp(disp_ol, 'false') == 1
end

% closed loop--------------------------------------------------------------
if strcmp(disp_cl, 'true') == 1

    figure
    subplot(5,1,1)
    hold on
    plot(t_sim, Vy, DisplayName='TruckSim')
    plot(t_sim, Vy_est_cl(1:end-1), DisplayName='Model')
    hold off
    title('Tractor Lateral Velocity $V_{y}$', 'Interpreter','latex')
    ylabel('m/s')
    legend(Location='best')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)
    
    subplot(5,1,2)
    hold on
    plot(t_sim, rad2deg(yaw_rate), DisplayName='TruckSim')
    plot(t_sim, rad2deg(yaw_rate_est_cl(1:end-1)), DisplayName='Model')
    hold off
    title('Yaw Rate $\dot{\psi}$', 'Interpreter', 'latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,3)
    hold on
    plot(t_sim, rad2deg(yaw), DisplayName='TruckSim')
    plot(t_sim, rad2deg(yaw_est_cl(1:end-1)), DisplayName='Model')
    hold off
    title('Yaw $\psi$', 'Interpreter','latex')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)
   
    subplot(5,1,4)
    hold on
    plot(t_sim, rad2deg(hitch_rate), DisplayName='hitch rate')
    plot(t_sim, rad2deg(hitch_rate_est_cl(1:end-1)), DisplayName='hitch rate')
    hold off
    title('Hitch Rate $\dot{\gamma}$', 'Interpreter','latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,5)
    hold on
    plot(t_sim, rad2deg(hitch), DisplayName='hitch')
    plot(t_sim, rad2deg(hitch_est_cl(1:end-1)), DisplayName='hitch')
    hold off
    title('Hitch $\gamma$', 'Interpreter','latex')
    xlabel('Time [s]')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

elseif strcmp(disp_cl, 'false') == 1
end

if strcmp(disp_cl, 'true') == 1
    
    % position solution comparison
    figure
    hold on
    plot(Xo, Yo, DisplayName='TruckSim')
    plot(pos_X_ol_est, pos_Y_ol_est, DisplayName='Model')
    plot(pos_X_cl_est, pos_Y_cl_est, DisplayName='KF')
    hold off
    title('Position Solutions')
    xlabel('Global X')
    ylabel('Global Y')
    legend(Location='best')
    grid
    set(gcf, 'color','w')

elseif strcmp(disp_cl, 'false') == 1
end