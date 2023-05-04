%% RLS Cornering Stiffness Estimator with Forgetting Factor 
function rls_cs = cornStiff(ts_data, Ay, Vy, x_derv, yaw_accel, hitch_accel, ...
            Vy_ol, yaw_rate_ol, yaw_ol, hitch_rate_ol,hitch_ol,  trunc_axle)
% Author: 
%           Tahn Thawainin, AU GAVLAB
%
% Description: 
%           A function to estimate tire cornering stiffnesses for a
%           5-axle tractor using a LS/RLS with forgetting factor
%
% Inputs: 
%           data - TruckSim data set (struct)
% Outputs: 
%           rls_cs - Cornering stiffness estimator data set (SI)

%% Load TruckSim Data

% simulation time
t_sim = ts_data.T_Event;

% average L1 and R1 steer angles(rad)
steer_ang = deg2rad((ts_data.Steer_L1 + ts_data.Steer_R1)/2);

% longitudinal velocity (m/s)
Vx = ts_data.Vx*(1e3/3600);

% lateral acceleration trailer
Ay2 = 9.81*ts_data.Ay_2;

% yaw (rad)
yaw = deg2rad(ts_data.Yaw);

% yaw rate (rad/s)
yaw_rate = deg2rad(ts_data.AVz);

% hitch (rad)
hitch = deg2rad(ts_data.Art_H);
% hitch = hotch_ol;

% hitch rated (rad/s)
hitch_rate = deg2rad(ts_data.ArtR_H);

% axle slip angles (rad)
sa1 = -deg2rad((ts_data.AlphaL1i + ts_data.AlphaR1i)./2);
sa2 = -deg2rad((ts_data.AlphaL2i + ts_data.AlphaR2i)./2);
sa3 = -deg2rad((ts_data.AlphaL3i + ts_data.AlphaR3i)./2);
sa4 = -deg2rad((ts_data.AlphaL4i + ts_data.AlphaR4i)./2);
sa5 = -deg2rad((ts_data.AlphaL5i + ts_data.AlphaR5i)./2);

% truncated slip angles (rad)
sa23 = (sa2 + sa3)./2;
sa45 = (sa4 + sa5)./2;

%% Vehicle Parameters
vp = vehParams();

%% Measurement from TruckSim

% measurment noise condition
% 0 - no added noise
% 1 - add measurement noise
meas_noise = 0;

if meas_noise == 0

    % lateral acceleration (m/s^2) - for cornering stiffness RLS
    Ay_meas = Ay;

elseif meas_noise == 1
    
    % lat accel measurement noise STD
    sigma_Ay = 0.1;
    
    for i = 1:length(t_sim)

    % lateral acceleration (m/s^2) - for cornering stiffness RLS
    Ay_meas(i) = (Ay(i) + Ay2(i)) + sigma_Ay*randn;

    end
end

%% Transients

% TS transients VS estimated transients
% 0 - No transients
% 1 - TS transients
% 2 - Estimated transients
trans_status = 2;

if trans_status == 1

yaw_accel = yaw_accel;
hitch_accel = hitch_accel;

elseif trans_status == 2

yaw_accel = x_derv(2,:);
hitch_accel = x_derv(4,:);
hitch = hitch_ol;

end

% trainsients
for i = 1:length(t_sim)
    t1(i) = -vp.m_t2*yaw_accel(i)*(vp.c + vp.d*cos(hitch(i)))/(vp.m_t1 + vp.m_t2);
    t2(i) = -vp.m_t2*hitch_accel(i)*vp.d*cos(hitch(i))/(vp.m_t1 + vp.m_t2);
end

%% Inputs

% Vy = Vy_ol;
% yaw_rate = yaw_rate_ol;
% yaw = yaw_ol;
% hitch_rate = hitch_rate_ol;
% hitch = hitch_ol;

%% RLS

% forgetting factor
lambda = 0.997;

% initialize
if trunc_axle == 0
x_init =  [3.2e5;
           1.5e5;
           1.5e5;
           1e5;
           1e5];

P_init = [1e5, 0, 0, 0, 0;...
          0, 1e2, 0, 0, 0;...
          0, 0, 1e2, 0, 0;...
          0, 0, 0, 1e2, 0;...
          0, 0, 0, 0, 1e2];

elseif trunc_axle == 1

x_init =  [3e5;...
           1.5e5;...
           1e5];

P_init = [1e1, 0, 0;...
          0, 1e1, 0;...
          0, 0, 1e1];...
end

x = x_init;
P = P_init;

% RLS
for i = 1:length(t_sim)
    
    if trunc_axle == 0

        if trans_status == 0

        % measurements
        y = (vp.m_t1 + vp.m_t2)*(Ay_meas(i));

        elseif trans_status == 1 || trans_status == 2

        % measurements
        y = (vp.m_t1 + vp.m_t2)*(Ay_meas(i) + t1(i) + t2(i));
        end

%         % scaling vector
        H = [steer_ang(i)*cos(steer_ang(i)) - cos(steer_ang(i))*((Vy(i) + vp.a*yaw_rate(i))/Vx(i)),...
             (-Vy(i) + vp.b1*yaw_rate(i))/Vx(i),...
             (-Vy(i) + vp.b2*yaw_rate(i))/Vx(i),...
             ((-cos(hitch(i))*Vy(i) + cos(hitch(i))*(vp.c + vp.f1*cos(hitch(i)))*yaw_rate(i) ...
                 + cos(hitch(i))^2*vp.f1*hitch_rate(i))/Vx(i)) + hitch(i)*cos(hitch(i)), ...
             ((-cos(hitch(i))*Vy(i) + cos(hitch(i))*(vp.c + vp.f2*cos(hitch(i)))*yaw_rate(i)...
                 + cos(hitch(i))^2*vp.f2*hitch_rate(i))/Vx(i)) + hitch(i)*cos(hitch(i))];
        
                % scaling vector
%         H = [steer_ang(i)*cos(steer_ang(i)) - cos(steer_ang(i))*((Vy(i) + vp.a*yaw_rate(i))/Vx(i)),...
%              (-Vy(i) + vp.b1*yaw_rate(i))/Vx(i),...
%              (-Vy(i) + vp.b2*yaw_rate(i))/Vx(i),...
%              ((-Vy(i) + (vp.c + vp.f1)*yaw_rate(i) + vp.f1*hitch_rate(i))/Vx(i)) + hitch(i), ...
%              ((-Vy(i) + (vp.c + vp.f2)*yaw_rate(i) + vp.f2*hitch_rate(i))/Vx(i)) + hitch(i)];

%         % TRUCKSIM SLIP ANGLES
%         H = [sa1(i), sa2(i), sa3(i), sa4(i), sa5(i)];
        
        % update gain matrix
        L = P*H'/(lambda + H*P*H');
        
        % update states
        innov = (y - H*x);
        x = x + L*(innov);
        
        % Update covariance matrix
        P = (eye(5) - L*H)*P*(1/lambda);
        
        %   direct least squares-------------------------------------------
        %     innov = 1;
        %     while abs(innov) > 0.1
        % 
        % %     update gain matrix
        %     L = P*H'./(lambda + H*P*H');
        % 
        % %     update states
        %     innov = (y - H*phi);
        %     x= x + pinv(H)*innov;
        
        %     end
        
        % siphon variables-------------------------------------------------
        
        % Cornering stiffnesses 
        rls_cs(i).C1 = x(1); 
        rls_cs(i).C2 = x(2);
        rls_cs(i).C3 = x(3);
        rls_cs(i).C4 = x(4);
        rls_cs(i).C5 = x(5);
        
        % slip angles
        rls_cs(i).sa1 = H(1);
        rls_cs(i).sa2 = H(2);
        rls_cs(i).sa3 = H(3);
        rls_cs(i).sa4 = H(4);
        rls_cs(i).sa5 = H(5);

        % gain
        rls_cs(i).L1 = L(1);
        rls_cs(i).L2 = L(2);
        rls_cs(i).L3 = L(3);
        rls_cs(i).L4 = L(4);
        rls_cs(i).L5 = L(5);

        % innovation
        rls_cs(i).innov = innov;

        % covariance
        rls_cs(i).P = P;

    elseif trunc_axle == 1

        % measurements
        y = (vp.m_t1 + vp.m_t2)*Ay_meas(i);

        % truncated distances
        b_trunc = vp.b1 + (vp.b2-vp.b1)/2;
        f_trunc = vp.f1 + (vp.f2-vp.f1)/2;
        
        % truncated observation matrix
        H = [steer_ang(i)*cos(steer_ang(i)) - cos(steer_ang(i))*((Vy(i)+vp.a*yaw_rate(i))/Vx(i)),...
             (-Vy(i) + b_trunc*yaw_rate(i))/Vx(i),...
             (-cos(hitch(i))*Vy(i) + cos(hitch(i))*(vp.c + f_trunc*cos(hitch(i)))*yaw_rate(i)...
             + cos(hitch(i))^2*f_trunc*hitch_rate(i))/Vx(i) + hitch(i)*cos(hitch(i))]; ...

        % TRUCKSIM SLIP ANGLES
%         H = [sa1(i), sa23(i), sa45(i)];

        % update gain matrix
        L = P*H'/(lambda + H*P*H');
        
        % update states
        innov = (y - H*x);
        x = x + L*(innov);
        
        % Update covariance matrix
        P = (eye(3) - L*H)*P*(1/lambda);
        
        %   direct least squares-------------------------------------------
        %     innov = 1;
        %     while abs(innov) > 0.1
        % 
        % %     update gain matrix
        %     L = P*H'./(lambda + H*P*H');
        % 
        % %     update states
        %     innov = (y - H*phi);
        %     x= x + pinv(H)*innov;
        %     
        % %     Update covariance matrix
        %     P = (eye(3) - L*H)*P*(1/lambda);
        
        %     end
        
        % siphon variables-------------------------------------------------

        % Cornering stiffnesses 
        rls_cs(i).C1 = x(1); 
        rls_cs(i).C23 = x(2);
        rls_cs(i).C45 = x(3);

        % slip angles
        rls_cs(i).sa1 = H(1);
        rls_cs(i).sa23 = H(2);
        rls_cs(i).sa45 = H(3);

        % gain
        rls_cs(i).L1 = L(1);
        rls_cs(i).L2 = L(2);
        rls_cs(i).L3 = L(3);

        % innovation
        rls_cs(i).innov = innov;

        % covariance
        rls_cs(i).P = P;
    end
end
end