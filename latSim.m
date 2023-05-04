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
ts_data = load('Run114_wideturn1.mat');

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
% steer_ang(1:100) = 0;

% steer ang test
steer_ang_test = deg2rad(5);

% longitudinal velocity (m/s)
Vx = ts_data.Vx*(1e3/3600);

% lateral velocity (m/s)
Vy = ts_data.Vy*(1e3/3600);

% lateral acceleration tractor
Ay = 9.81*ts_data.Ay;

% yaw (rad)
yaw = deg2rad(ts_data.Yaw);

% yaw rate (rad/s)
yaw_rate = deg2rad(ts_data.AVz);

% derive yaw accel (rad/s^2)
yaw_accel(1) = 0;
for i = 2:length(t_sim)
    yaw_accel(i) = (yaw_rate(i) - yaw_rate(i-1))/dt;
end

% hitch (rad)
hitch = deg2rad(ts_data.Art_H);

% hitch rated (rad/s)
hitch_rate = deg2rad(ts_data.ArtR_H);

% derive hitch accel (rad/s^2) 
hitch_accel(1) = 0;
for i = 2:length(t_sim)
    hitch_accel(i) = (hitch_rate(i) - hitch_rate(i-1))/dt;
end

% transients
for i = 1:length(t_sim)
    t1(i) = (-vp.m_t2*yaw_accel(i)*(vp.c + vp.d*cos(hitch(i))))/(vp.m_t1 + vp.m_t2);
    t2(i) = (-vp.m_t2*hitch_accel(i)*vp.d*cos(hitch(i)))/(vp.m_t1 + vp.m_t2);
end

% total acceleration
Ay_tot = Ay + t1' + t2';

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

%% Open-Loop Dynamic Model Simulation

% initialize states
Vy_ol = 0;
yaw_rate_ol = 0;
yaw_ol = 0;
hitch_rate_ol = 0;
hitch_ol = 0;

% initialize position
pos_X_ol = 0;
pos_Y_ol = 0;

x_init = [Vy_ol;
          yaw_rate_ol;
          yaw_ol;
          hitch_rate_ol;
          hitch_ol];

x = x_init;

% tire stiffness values
% CS = [ltm.A1_cs, ltm.A23_cs, ltm.A23_cs, ltm.A45_cs, ltm.A45_cs]; % truth?
CS = [4.e5, 2.5e5, 2.5e5, 1.5e5, 1.5e5];
% CS = [2.5e5, 2.5e5, 2.5e5, 1e5, 1e5];

Vx_const = 25;

for i = 1:length(t_sim)

% input
u = steer_ang(i);

% update dynamics
lat_ol(i) = latModel(steer_ang(i), Vx(i), dt, CS);

% simulate dynamics
xd = lat_ol(i).sysc.A*x + lat_ol(i).sysc.B*u;

% update states
x = x + xd.*dt;

% position integration
pos_X_ol = pos_X_ol + (Vx(i)*cos(x(3)) - x(1)*sin(x(3)))*dt;
pos_Y_ol = pos_Y_ol + (Vx(i)*sin(x(3)) + x(1)*cos(x(3)))*dt;

% siphon variables---------------------------------------------------------

% state derivatives
x_derv(:,i) = xd;

% state estimates
Vy_ol(i,:) = x(1);
yaw_rate_ol(i,:) = x(2);
yaw_ol(i,:) = x(3);
hitch_rate_ol(i,:) = x(4);
hitch_ol(i,:) = x(5);

% position 
pos_X_ol_est(i,:) = pos_X_ol;
pos_Y_ol_est(i, :) = pos_Y_ol;

end

Ay_ol = x_derv(1,:);

%% Pre-Process Tire Stiffness Estimattion
% Uses clean signals from TruckSim

% truncate axle condition
% 0 - DONT truncate axles
% 1 - truncate axles
trunc_axle = 0;

% call estimator
rls_cs = cornStiff(ts_data, Ay, Vy, x_derv, yaw_accel, hitch_accel, ...
            Vy_ol, yaw_rate_ol, yaw_ol, hitch_rate_ol, hitch_ol, trunc_axle);

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
C23_est = extractfield(rls_cs, 'C23');
C45_est = extractfield(rls_cs, 'C45');

sa1_est = extractfield(rls_cs, 'sa1');
sa23_est = extractfield(rls_cs, 'sa23');
sa45_est = extractfield(rls_cs, 'sa45');

end

%% Measurements from TruckSim

% measurment noise condition
% 0 - no added noise
% 1 - add measurement noise
meas_noise = 0;

if meas_noise == 0

    % lateral acceleration (m/s^2)
    Ay_meas = Ay;    
    % yaw rate
    yaw_rate_meas = yaw_rate;

elseif meas_noise == 1
    
    % lat accel measurement noise STD
    sigma_Ay = 0.05;
    % yaw rate measurement noise STD
    sigma_yr = 0.001;
    % hitch rate measurement noise STD
    sigma_hr = 0.001;
    
    % noise vector
    n_Ay = sigma_Ay.*randn(length(t_sim),1);
    n_yr = sigma_yr.*randn(length(t_sim),1);
    n_hr = sigma_hr.*randn(length(t_sim),1);
    
    % lateral acceleration (m/s^2) 
    Ay_meas = Ay + sigma_Ay.*randn(length(t_sim),1);
    % yaw rate
    yaw_rate_meas = yaw_rate + sigma_yr.*randn(length(t_sim),1);
    % hitch rate
    hitch_rate_meas = hitch_rate + sigma_hr.*randn(length(t_sim),1);

    % signal to noise ratio
    SNR_Ay = snr(Ay, n_Ay);
    SNR_yr = snr(yaw_rate, n_yr);
    SNR_hr = snr(hitch_rate, n_hr);

end


%% Closed Loop Kalman Filter

% PROCESS AND MEASUREMENT NOISE--------------------------------------------

Q =  [10, 0, 0, 0, 0;
      0, 10, 0, 0, 0;
      0, 0, 10, 0, 0;
      0, 0, 0, 10, 0;
      0, 0, 0, 0, 10];

% INITIALIZE---------------------------------------------------------------

% initialize states
Vy_cl = 0;
yaw_rate_cl = 0;
yaw_cl = 0;
hitch_rate_cl = 0;
hitch_cl = 0;

% initialize position
pos_X_cl = 0;
pos_Y_cl = 0;

x_init = [Vy_cl;
          yaw_rate_cl;
          yaw_cl;
          hitch_rate_cl;
          hitch_cl];

P_init = [10, 0, 0, 0, 0;...
          0, 10, 0, 0, 0;...
          0, 0, 10, 0, 0;...
          0, 0, 0, 10, 0;...
          0, 0, 0, 0, 10];

x = x_init;
P = P_init;

% TIME UPDATE--------------------------------------------------------------

for i = 1:length(t_sim)

% input
u = steer_ang(i);

% update dynamics
lat_ol(i) = latModel(steer_ang(i), Vx(i), dt, CS);

% simulate dynamics
x = lat_ol(i).sysd.A*x + lat_ol(i).sysd.B*u;

% priori covariance
P = lat_ol(i).sysd.A*P*lat_ol(i).sysd.A' + Q;

% MEASUREMENT UPDATE-------------------------------------------------------
% use hitch rate meas condition
% 0 - DONT use
% 1 - USE
h_meas = 0;

% no puesdo measurment
if h_meas == 0 

    % measurement noise
    if meas_noise == 0
    R = [0.01, 0;
         0, 0.01];
    elseif meas_noise == 1
    R = [1.5, 0;
         0, 1.5];
    end

% measurements
y = [Ay_meas(i);
     yaw_rate_meas(i)];

% observation matrix
H = [0, Vx(i), 0, 0, 0;
     0, 1, 0, 0, 0];

elseif h_meas == 1

    % measurement noise
    if meas_noise == 0
    R = [0.01, 0, 0;
         0, 0.01, 0;
         0, 0, 0.01];
    elseif meas_noise == 1
    R = [1, 0, 0;
         0, 1, 0;
         0, 0, 3];
    end

% trailer split axle length
L = (vp.f1 + vp.f2)/2;

% measurements
y = [Ay_meas(i);
     yaw_rate_meas(i);
     hitch_rate_meas(i)];

% observation matrix with puesdo measurement
H = [0, Vx(i), 0, 0, 0;
     0, 1, 0, 0, 0;
     0, 0, 0, 1, 0];
end

% kalman gain
K = P*H'/(H*P*H' + R);

% posteriori covariance
P = (eye(5) - K*H)*P;

% measurement prediction innovation
innov = (y - H*x);

% state correction update
x = x + K*(innov);

% nonlienar transform------------------------------------------------------
d_v = -sin(x(3))*dt;
d_y = (-Vx(i)*sin(x(3)) - x(1)*cos(x(3)))*dt;

% nonlinear observation
H_lin = [d_v, 0, 0, 0, 0;
         0, 0, d_y, 0, 0];

% transformed covariance
P_trans = H_lin*P*H_lin';

% position integration
pos_X_cl = pos_X_cl + (Vx(i)*cos(x(3)) - x(1)*sin(x(3)))*dt;
pos_Y_cl = pos_Y_cl + (Vx(i)*sin(x(3)) + x(1)*cos(x(3)))*dt;

% SIPHON VARIABLES---------------------------------------------------------

% state estimates
Vy_cl(i,:) = x(1);
yaw_rate_cl(i,:) = x(2);
yaw_cl(i,:) = x(3);
hitch_rate_cl(i,:) = x(4);
hitch_cl(i,:) = x(5);

% position estimates
pos_X_cl_est(i,:) = pos_X_cl;
pos_Y_cl_est(i,:) = pos_Y_cl;

% transformed covariance
covar_trans(i,:,:) = P_trans;

end

%% Error Analysis

% open loop error
X_err_ol = pos_X_ol_est - Xo;
Y_err_ol = pos_Y_ol_est - Yo;

% closed loop error
X_err_cl = pos_X_cl_est - Xo;
Y_err_cl = pos_Y_cl_est - Yo;

% norms
for i = 1:length(t_sim)
    norm_ts(i) = norm(Xo(i), Yo(i));
    norm_ol(i) = norm(pos_X_ol_est(i), pos_Y_ol_est(i));
    norm_cl(i) = norm(pos_X_cl_est(i), pos_Y_cl_est(i));
end

%% Interface

% colors
green = [0.4660 0.6740 0.1880];
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];

% display tire model
disp_tire_model = 'false';
% display tire pre-process tire siffness estimator
disp_pre_ts = 'false';
% display true states
disp_true_states = 'false';
% dispaly open loop dynamics
disp_ol = 'false';
% display closed loop dynamics
disp_cl = 'false';
% disply position solution
disp_pos = 'true';
% display comparision plots
disp_comp = 'true';
% display error plots
disp_err = 'true';

% tire model---------------------------------------------------------------
if strcmp(disp_tire_model, 'true') == 1
    
    % overlap figure
    % 0 - DONT overlap
    % 1 - overlap
    fig_overlap = 0;

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
    plot(t_sim, C5_est, DisplayName='C5')
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
    plot(t_sim, C1_est(end).*sa1_est, DisplayName='RLS')
%     plot(t_sim, ltm.A1_cs*sa1_est, DisplayName='RLS')
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
    plot(t_sim, C2_est(end).*sa2_est, DisplayName='RLS')
%     plot(t_sim, ltm.A23_cs*sa2_est, DisplayName='RLS')
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
    plot(t_sim, C3_est(end).*sa3_est, DisplayName='RLS')
%     plot(t_sim, ltm.A23_cs*sa3_est, DisplayName='RLS')
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
    plot(t_sim, C4_est(end).*sa4_est, DisplayName='RLS')
%     plot(t_sim, ltm.A45_cs*sa4_est, DisplayName='RLS')
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
    plot(t_sim, C5_est(end).*sa5_est, DisplayName='RLS')
%     plot(t_sim, ltm.A45_cs*sa5_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 5')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')

    % --- truncated axles ------
    elseif trunc_axle == 1

    % stiffness estimation
    figure
    hold on
    plot(t_sim, C1_est, DisplayName='C1')
    plot(t_sim, C23_est, DisplayName='C23')
    plot(t_sim, C45_est, DisplayName='C45')
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
    plot(t_sim, C1_est(end).*sa1_est, DisplayName='RLS')
%     plot(t_sim, ltm.A1_cs*sa1_est, DisplayName='RLS')
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
    plot(t_sim, ltm.A23_cs*sa23, DisplayName='Sim Truth')
    plot(t_sim, C23_est(end).*sa23_est, DisplayName='RLS')
%     plot(t_sim, ltm.A23_cs*sa2_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 2')
    xlabel('Time [s]')
    ylabel('N')
    legend(Location='best')
    grid
    set(gcf, 'color','w')

    % Lateral Force comparison
    figure
    hold on
    plot(t_sim, ts_data.Fy_A4, DisplayName='TS')
    plot(t_sim, ltm.A45_cs*sa45, DisplayName='Sim Truth')
    plot(t_sim, C45_est(end).*sa45_est, DisplayName='RLS')
%     plot(t_sim, ltm.A45_cs*sa4_est, DisplayName='RLS')
    hold off
    title('Fy at Axle 4')
    xlabel('Time [s]')
    ylabel('N')
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
    plot(t_sim, Vy_ol, DisplayName='Model')
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
    plot(t_sim, rad2deg(yaw_rate_ol), DisplayName='Model')
    hold off
    title('Yaw Rate $\dot{\psi}$', 'Interpreter', 'latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,3)
    hold on
    plot(t_sim, rad2deg(yaw), DisplayName='TruckSim')
    plot(t_sim, rad2deg(yaw_ol), DisplayName='Model')
    hold off
    title('Yaw $\psi$', 'Interpreter','latex')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)
   
    subplot(5,1,4)
    hold on
    plot(t_sim, rad2deg(hitch_rate), DisplayName='hitch rate')
    plot(t_sim, rad2deg(hitch_rate_ol), DisplayName='hitch rate')
    hold off
    title('Hitch Rate $\dot{\gamma}$', 'Interpreter','latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,5)
    hold on
    plot(t_sim, rad2deg(hitch), DisplayName='hitch')
    plot(t_sim, rad2deg(hitch_ol), DisplayName='hitch')
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
    
    % closed loop KF plots
    figure
    subplot(5,1,1)
    hold on
    plot(t_sim, Vy, DisplayName='TruckSim')
    plot(t_sim, Vy_cl, DisplayName='KF')
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
    plot(t_sim, rad2deg(yaw_rate_cl), DisplayName='KF')
    hold off
    title('Yaw Rate $\dot{\psi}$', 'Interpreter', 'latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,3)
    hold on
    plot(t_sim, rad2deg(yaw), DisplayName='TruckSim')
    plot(t_sim, rad2deg(yaw_cl), DisplayName='KF')
    hold off
    title('Yaw $\psi$', 'Interpreter','latex')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)
   
    subplot(5,1,4)
    hold on
    plot(t_sim, rad2deg(hitch_rate), DisplayName='hitch rate')
    plot(t_sim, rad2deg(hitch_rate_cl), DisplayName='hitch rate')
    hold off
    title('Hitch Rate $\dot{\gamma}$', 'Interpreter','latex')
    ylabel('deg/s')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

    subplot(5,1,5)
    hold on
    plot(t_sim, rad2deg(hitch), DisplayName='hitch')
    plot(t_sim, rad2deg(hitch_cl), DisplayName='hitch')
    hold off
    title('Hitch $\gamma$', 'Interpreter','latex')
    xlabel('Time [s]')
    ylabel('deg')
    grid
    set(gcf, 'color','w')
%     set(gca, 'fontsize', 16)

elseif strcmp(disp_cl, 'false') == 1
end

% position solutions-------------------------------------------------------
if strcmp(disp_pos, 'true') == 1
    
    % position solution comparison
    figure
    hold on
    plot(Xo, Yo, color = green, DisplayName='TruckSim', color = green, LineWidth=1.5)
    plot(pos_X_ol_est, pos_Y_ol_est, DisplayName='Model', color = blue, LineWidth=1.5)
    plot(pos_X_cl_est, pos_Y_cl_est, DisplayName='KF', color = orange, LineWidth=1.5)
    hold off
    title('Position Solution')
    xlabel('Global X [m]')
    ylabel('Global Y [m]')
    leg = legend;
    set(leg,'Interpreter','latex');
    grid
    set(gcf, 'color','w')

elseif strcmp(disp_pos, 'false') == 1
end

% comparision plots--------------------------------------------------------
if strcmp(disp_comp, 'true') == 1
    
    % lateral velocity
    figure
    hold on
    plot(t_sim, Vy, color = green, LineWidth=1.5)
    plot(t_sim, Vy_ol, color = blue, LineWidth=1.5)
    plot(t_sim, Vy_cl, color = orange, LineWidth=1.5)
    hold off
    title('Lateral Velocity $V_{y}$', 'Interpreter','latex')
    xlabel('Time [s]', 'Interpreter','latex')
    ylabel('$m/s$', 'Interpreter','latex')
    legend('TruckSim','Model', 'KF', 'Interpreter', 'latex')
    grid
    set(gcf, 'color', 'w')

    % yaw
    figure
    hold on
    plot(t_sim, rad2deg(yaw), color = green, DisplayName='TruckSim', LineWidth=1.5)
    plot(t_sim, rad2deg(yaw_ol), DisplayName='Model', color = blue, LineWidth=1.5)
    plot(t_sim, rad2deg(yaw_cl), DisplayName='KF', color = orange, LineWidth=1.5)
    hold off
    title('Yaw $\psi$', 'Interpreter','latex')
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$Deg$','Interpreter','latex')
    leg = legend;
    set(leg,'Interpreter','latex');
    grid
    set(gcf, 'color', 'w')

    % yaw rate
    figure
    hold on
    plot(t_sim, rad2deg(yaw_rate), color = green, DisplayName='TruckSim', LineWidth=1.5)
    plot(t_sim, rad2deg(yaw_rate_ol), DisplayName='Model', color = blue, LineWidth=1.5)
    plot(t_sim, rad2deg(yaw_rate_cl), DisplayName='KF', color = orange, LineWidth=1.5)
    hold off
    title('Yaw Rate $\dot{\psi}$', 'Interpreter','latex')
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$Deg/s$','Interpreter','latex')
    leg = legend;
    set(leg,'Interpreter','latex');
    grid
    set(gcf, 'color', 'w')

    % hitch
    figure
    hold on
    plot(t_sim, rad2deg(hitch), color = green, DisplayName='TruckSim', LineWidth=1.5)
    plot(t_sim, rad2deg(hitch_ol), DisplayName='Model', color = blue, LineWidth=1.5)
    plot(t_sim, rad2deg(hitch_cl), DisplayName='KF', color = orange, LineWidth=1.5)
    hold off
    title('Hitch $\gamma$', 'Interpreter','latex')
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$Deg$','Interpreter','latex')
    leg = legend;
    set(leg,'Interpreter','latex');
    grid
    set(gcf, 'color', 'w')

    % hitch rate
    figure
    hold on
    plot(t_sim, rad2deg(hitch_rate), color = green, DisplayName='TruckSim', LineWidth=1.5)
    plot(t_sim, rad2deg(hitch_rate_ol), DisplayName='Model', color = blue, LineWidth=1.5)
    plot(t_sim, rad2deg(hitch_rate_cl), DisplayName='KF', color = orange, LineWidth=1.5)
    hold off
    title('Hitch Rate $\dot{\gamma}$', 'Interpreter','latex')
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$Deg/s$','Interpreter','latex')
    leg = legend;
    set(leg,'Interpreter','latex');
    grid
    set(gcf, 'color', 'w')

    % lateral acceleration
    figure
    hold on
    plot(t_sim, Ay, color = green, LineWidth=1.5)
    plot(t_sim, yaw_rate_ol.*Vx, color = blue, LineWidth=1.5)
    plot(t_sim, yaw_rate_cl.*Vx, color = orange, LineWidth=1.5)
    hold off
    legend('Model $\dot{\psi}V_{x}$', 'KF $\dot{\psi}V_{x}$', 'TruckSim',...
           'Interpreter', 'latex')
    title('Lateral Acceleration $A_{y}$', 'Interpreter','latex')
    xlabel('Time [s]','Interpreter','latex')
    ylabel('$m/s^{2}$','Interpreter','latex')
    grid
    set(gcf, 'color', 'w')

elseif strcmp(disp_comp, 'false') == 1
end

% error plots--------------------------------------------------------------
if strcmp(disp_err, 'true') == 1
    
    figure
    hold on
    plot(t_sim, Y_err_ol, color = blue, LineWidth=1.5)
    plot(t_sim, Y_err_cl, color = orange, LineWidth=1.5)
    hold off
    title('Position Error')
    xlabel('Time [s]')
    ylabel('[m]')
    legend('Model','KF', 'Interpreter','latex')
    grid
    set(gcf, 'color', 'w')



elseif strcmp(disp_err, 'false') == 1
end