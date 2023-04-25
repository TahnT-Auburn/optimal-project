%% Lateral Tire Model
function ltm = latTireModel(ts_tiremodel, vert_load)

% Author: Tahn Thawainin, AU GAVLAB
%
% Description: A function to store the TruckSim lateral tire model. Includes
%              interpolated tire curves for specified vertical loads
%             
% Input(s): ts_tiremodel - (.csv file) internal TruckSim tire model
%           vert_load - array of vertical forces at each tire (T(1)-T(k))
%
% Output(s): ltm - lateral tire model data set

%% Load TruckSim Tire Model

tm = readtable(ts_tiremodel);

% convert table to array
ltm.tm_array = table2array(tm);

% variables
ltm.slip_ang = deg2rad(ltm.tm_array(2:end,1));     % slip angle (rad)
ltm.Fz = ltm.tm_array(1,2:end);           % vertical force
ltm.Fy1 = ltm.tm_array(2:end,2);          % lateral force at Fz1
ltm.Fy2 = ltm.tm_array(2:end,3);          % lateral force at Fz2
ltm.Fy3 = ltm.tm_array(2:end,4);          % lateral force at Fz3
ltm.Fy4 = ltm.tm_array(2:end,5);          % lateral force at Fz4
ltm.Fy5 = ltm.tm_array(2:end,6);          % lateral force at Fz5

% Fz vs Fy
Fz_vs_Fy = ltm.tm_array(1:end,2:end);

%% Tire Curve Interpolation

    for i = 1:length(ltm.slip_ang)
        
        % 1st tire curve---------------------------------------------------

        % lower limit
        low_lim = find(ltm.Fz < vert_load(1));
        low_lim = low_lim(end);

        % upper limit
        up_lim = find(ltm.Fz > vert_load(1));
        up_lim = up_lim(1);
        
        Fy_T1(i) = Fz_vs_Fy(i+1,low_lim) + (vert_load(1) - ltm.Fz(low_lim))*...
                   ((Fz_vs_Fy(i+1,up_lim) - Fz_vs_Fy(i+1,low_lim))/...
                   (ltm.Fz(up_lim) - ltm.Fz(low_lim)));

        % 2nd and 3rd tire cure--------------------------------------------

        % lower limit
        low_lim = find(ltm.Fz < vert_load(2));
        low_lim = low_lim(end);

        % upper limit
        up_lim = find(ltm.Fz > vert_load(2));
        up_lim = up_lim(1);
        
        Fy_T2(i) = Fz_vs_Fy(i+1,low_lim) + (vert_load(2) - ltm.Fz(low_lim))*...
                   ((Fz_vs_Fy(i+1,up_lim) - Fz_vs_Fy(i+1,low_lim))/...
                   (ltm.Fz(up_lim) - ltm.Fz(low_lim)));

        % 4th and 5th tire curve-------------------------------------------

        % lower limit
        low_lim = find(ltm.Fz < vert_load(3));
        low_lim = low_lim(end);

        % upper limit
        up_lim = find(ltm.Fz > vert_load(3));
        up_lim = up_lim(1);
        
        Fy_T3(i) = Fz_vs_Fy(i+1,low_lim) + (vert_load(3) - ltm.Fz(low_lim))*...
                   ((Fz_vs_Fy(i+1,up_lim) - Fz_vs_Fy(i+1,low_lim))/...
                   (ltm.Fz(up_lim) - ltm.Fz(low_lim)));
    end

    % siphone variables
    ltm.Fy_T1 = Fy_T1;
    ltm.Fy_T2 = Fy_T2;
    ltm.Fy_T3 = Fy_T3;
    
    % tire cornering stiffness
    ltm.tire1_cs = (ltm.Fy_T1(4) - ltm.Fy_T1(1))/...
                   (ltm.slip_ang(4) - ltm.slip_ang(1));

    ltm.tire2_cs = (ltm.Fy_T2(4) - ltm.Fy_T2(1))/...
                   (ltm.slip_ang(4) - ltm.slip_ang(1));

    ltm.tire3_cs = (ltm.Fy_T3(4) - ltm.Fy_T3(1))/...
                   (ltm.slip_ang(4) - ltm.slip_ang(1));
    
    % effective axle cornering stiffness
    ltm.A1_cs = 2*ltm.tire1_cs;
    ltm.A23_cs = 2*ltm.tire2_cs;
    ltm.A45_cs = 2*ltm.tire3_cs;

end