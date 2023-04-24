%% Longitudinal Velocity Controllers
function contr = lonContr(contrType)
% Author: Tahn Thawainin, AU GAVLAB
%
% Description: A function generate a P, PD, PI, or PID longitudinal velocity
%              controller gains
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
%
% Ouputs: contr - longitudinal controller data set (SI)

% Vehicle Parameters
vp = vehParams();

% Simulation proedure
pro = lonProc();

%% Controller Specs

% Proportional (P)---------------------------------------------------------

% time constant
contr.t_P = 5;

% Controller Gains

% Proportional (P)---------------------------------------------------------
if contrType == 0

% proportional gain
contr.Kp = (pro.M_eff/contr.t_P - pro.B_eff)*(1/pro.scale_factor);

% DC gain
contr.DC = (pro.scale_factor*contr.Kp)/(pro.B_eff + pro.scale_factor*contr.Kp);
else 
end

end