%% Nonlinear Recursive Least Squares Mass Estimation
%
% Description: Script to perform vehicle mass estimation using a nonlinear
%              recursive least squares algorithm
%
% Author: Tahn Thawainin, AU GAVLAB

clc
clear variables
close all

%% Load Data 

% simulation data set
data = load("u5a_p_sg_ct1.mat");

% vehicle parameters
vp = data.vp;
% procedure
pro = data.pro;
% simulation
sim = data.sim;

% Extract Simulation Variables
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

%% Measurements

% measurement variance
sigma_v = 0.01;

% accelerometer bias
b_a = 0;

% measurement
y = b_a + accel + sigma_v*randn(1,length(pro.t_sim));

%% Nonlinear Recursive Least Squares

% intialize
x_init = 1/1000;
P_init = 1e-6;

x = x_init;
P = P_init;

% measurement noise
R = 0.1;

for k = 1:length(pro.t_sim)
    
    % linearized observation matrix
    H = pro.scale_factor*T_eng(k) - pro.B_eff*vel(k)...
        - vp.u_rr*pro.g*cos(pro.grade(k)) - pro.g*sin(pro.grade(k)) ...
        - 0.5*pro.p*vp.cd*vp.front_area*vel(k)^2;

    % nonlinear observation matrix
    h = pro.scale_factor*T_eng(k) - pro.B_eff*vel(k)...
        - ((1/x(1)) - pro.M_i)*vp.u_rr*pro.g*cos(pro.grade(k)) ...
        - ((1/x(1)) - pro.M_i)*pro.g*sin(pro.grade(k)) ...
        - 0.5*pro.p*vp.cd*vp.front_area*vel(k)^2;
    
%     % measurement noise
%     if abs(y(k)) > 0.05
%         R = 1e-5;
%     elseif abs(y(k)) < 0.05
%         R = 1;
%     end
    
    % gain
    L = P*H'/(H*P*H' + R);

    % covariance
    P = (eye(1) - L*H)*P;

    % state estimate
    x = x + L*(y(k) - h*x);
    
    % linearized transform
    H_lin = -1/x^2;

    % transformed covariance
    P_trans = H_lin*P*H_lin';

    % siphon variables-----------------------------------------------------

    % effective mass estimate
    M_eff_est(k) = (1/x(1));
    % vehicle mass estimate
    M_veh_est(k) = (1/x(1)) - pro.M_i;
    % accelerometer bias estimate
%     b_a_est(k) = x(2);
    
    % observability 
    obsv(k,:) = H;
    % kalman gain
    gain(k,:) = L;
    % covariance matrix
    covar(k,:) = P;
    % transformed covariance matrix
    covar_trans(k,:) = P_trans;
end

    % transformed theoretical STD (P)
    theo_std_trans = sqrt(covar);
%% Observability Analysis

if rank(obsv) == length(x)
    obsv_status = 'Observable';
else
    obsv_status = ' Unobservable';
end

%% Monte Carlo Analysis
% monte carlo conidition
% 0 - DONT run monte carlo
% 1 - Run monte carlo
monte_carlo = 1;

% monte carlo length
mc_length = 100;

% run monte carlo
if monte_carlo == 1

    for j = 1:mc_length
        
        % measurement
        y = b_a + accel + sigma_v*randn(1,length(pro.t_sim));
        
        % Nonlinear Recursive Least Squares
        
        % intialize
        x_init = 1/1000;
        P_init = 1e-6;

        
        x = x_init;
        P = P_init;
        
        % measurement noise
        R = 0.1;
        
        for k = 1:length(pro.t_sim)
            
            % linearized observation matrix
            H = pro.scale_factor*T_eng(k) - pro.B_eff*vel(k)...
                - vp.u_rr*pro.g*cos(pro.grade(k)) - pro.g*sin(pro.grade(k)) ...
                - 0.5*pro.p*vp.cd*vp.front_area*vel(k)^2;
        
            % nonlinear observation matrix
            h = pro.scale_factor*T_eng(k) - pro.B_eff*vel(k)...
                - ((1/x(1)) - pro.M_i)*vp.u_rr*pro.g*cos(pro.grade(k)) ...
                - ((1/x(1)) - pro.M_i)*pro.g*sin(pro.grade(k)) ...
                - 0.5*pro.p*vp.cd*vp.front_area*vel(k)^2;

        %     % measurement noise
        %     if abs(y(k)) > 0.05
        %         R = 1e-5;
        %     elseif abs(y(k)) < 0.05
        %         R = 1;
        %     end
            
            % gain
            L = P*H'/(H*P*H' + R);
        
            % covariance
            P = (eye(1) - L*H)*P;
        
            % state estimate
            x = x + L*(y(k) - h*x);
        
            % siphon variables---------------------------------------------
            
            % estimated state (1/Meff)
            state(k) = x(1);
            % effective mass estimate
            M_eff_est(k) = (1/x(1));
            % vehicle mass estimate
            M_veh_est(k) = (1/x(1)) - pro.M_i;
            % accelerometer bias estimate
%             b_a_est(k) = x(2);
            
            % observability 
            obsv(k,:) = H;
            % kalman gain
            gain(k,:) = L;
            % covariance matrix
            covar(k,:) = P;
        end

        % monte carlo results----------------------------------------------
        
        % estimated state (1/Meff)
        state_mc(:,j) = state;

        % estimated state errors
        state_err_mc(:,j) = state - (1/pro.M_eff);
        
        % vehicle mass estimates
        m_veh_mc(:,j) = M_veh_est;

        % vehicle mass errors
        m_veh_err_mc(:,j) = M_veh_est - vp.m_veh;
    end
        
    % monte carlo analysis (mean and STD)
    state_mc_mean = mean(state_mc,2);
    state_err_mc_mean = mean(state_err_mc,2);
    m_mc_mean = mean(m_veh_mc,2);
    m_err_mc_mean = mean(m_veh_err_mc,2);
    
    state_mc_std = std(state_mc, 0,2);
    state_err_mc_std = std(state_mc, 0,2);
    m_mc_std = std(m_veh_mc,0,2);
    m_err_mc_std = std(m_veh_err_mc,0,2);
    
    % theoretical STD (P)
    theo_std = sqrt(covar(1,1,:));
    theo_std = (reshape(theo_std, 1,length(theo_std)));
    
elseif monte_carlo == 0
end
      
%% Interface

% colors
grayColor = [.7 .7 .7];

% display vehicle info
veh_info = 'false';
% display procedure info
pro_info = 'false';
% display simulation info
sim_info = 'false';
% display estimator info
est_info = 'true';
% display monte carlo info
mc_info = 'true';

% vehicle info-------------------------------------------------------------
if strcmp(veh_info, 'true') == 1
    
    disp('Vehicle Configuration:')

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

% procedure info-----------------------------------------------------------
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

% simulation information---------------------------------------------------
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
    plot(pro.t_sim, accel)
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

% estimator information----------------------------------------------------
if strcmp(est_info, 'true') == 1
    
    % estimator specs
    disp('Estimator Specs:')
    
    disp('Observability:')
    disp(obsv_status)

    disp('Measurement Noise:')
    disp(R)

    % NRLS comparison plot
    figure
    set(gcf,'color','w')
    hold on
    plot(pro.t_sim, vp.m_veh*ones(1,length(pro.t_sim)), '--', LineWidth = 1.5)
    plot(pro.t_sim, M_veh_est, LineWidth = 1.5)
    hold off
    title('Vehicle Mass Estimate')
    xlabel('Time [s]')
    ylabel('Mass [Kg]')
    ylim([(vp.m_veh - 0.03*vp.m_veh),(vp.m_veh + 0.03*vp.m_veh)])
    grid
    legend('Truth', 'NRLS', Location='best')
    
%     % estimated bias plot
%     figure
%     set(gcf, 'Color', 'w')
%     hold on
%     plot(pro.t_sim, b_a*ones(1,length(pro.t_sim)), '--', LineWidth = 1.5)
%     plot(pro.t_sim, b_a_est, LineWidth = 1.5)
%     hold off
%     title('Accelerometer Bias Estimate')
%     xlabel('Time [s]')
%     ylabel('Bias [m/s]')
%     grid
%     legend('Truth', 'NRLS', Location='best')

elseif strcmp(est_info, 'false') == 1
end

% vehicle info-------------------------------------------------------------
if strcmp(mc_info, 'true') == 1

    figure
    subplot(2,1,1)
    hold on
    plot(pro.t_sim, m_mc_mean + 3.*m_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, m_mc_mean + 3.*theo_std_trans, '--k', LineWidth=2)
    plot(pro.t_sim, m_mc_mean - 3.*m_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, m_mc_mean - 3.*theo_std_trans, '--k', LineWidth=2)
%     plot(pro.t_sim, m_veh_mc ,'color', grayColor)
    hold off
    title('Vehicle Mass Estimates')
    ylabel('Estimates [Kg]')
    ylim([(vp.m_veh - 0.015*vp.m_veh),(vp.m_veh + 0.015*vp.m_veh)])
    legend('Experimental STD') % 'Theoretical STD')
    subplot(2,1,2)
    hold on
    plot(pro.t_sim, m_err_mc_mean + 3.*m_err_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, m_err_mc_mean + 3.*theo_std_trans, '--k', LineWidth=2)
    plot(pro.t_sim, m_err_mc_mean - 3.*m_err_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, m_err_mc_mean - 3.*theo_std_trans, '--k', LineWidth=2)
%     plot(pro.t_sim, m_veh_err_mc ,'color', grayColor)
    hold off
    title('Vehicle Mass Estimate Errors')
    ylabel('Estimate Errors [kg]')   
    xlabel('Time [s]')
    ylim([-200,200])
    set(gcf, 'color', 'w')
    
    figure
    subplot(2,1,1)
    hold on
    plot(pro.t_sim, state_mc_mean + 3.*state_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, state_mc_mean + 3.*theo_std', '--k', LineWidth=2)
    plot(pro.t_sim, state_mc_mean - 3.*state_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, state_mc_mean - 3.*theo_std', '--k', LineWidth=2)
    plot(pro.t_sim, state_mc ,'color', grayColor)
    hold off
    title('State Estimates')
    ylabel('Estimates [1/Kg]')
%     ylim([(vp.m_veh - 0.015*vp.m_veh),(vp.m_veh + 0.015*vp.m_veh)])
    legend('Experimental STD', 'Theoretical STD')
    subplot(2,1,2)
    hold on
    plot(pro.t_sim, state_err_mc_mean + 3.*state_err_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, state_err_mc_mean + 3.*theo_std', '--k', LineWidth=2)
    plot(pro.t_sim, state_err_mc_mean - 3.*state_err_mc_std, 'r', LineWidth=2)
    plot(pro.t_sim, state_err_mc_mean - 3.*theo_std', '--k', LineWidth=2)
    plot(pro.t_sim, state_err_mc ,'color', grayColor)
    hold off
    title('State Estimate Errors')
    ylabel('Estimate Errors [1/kg]')   
    xlabel('Time [s]')
%     ylim([-200,200])
    set(gcf, 'color', 'w')
elseif strcmp(mc_info, 'false') == 0
end