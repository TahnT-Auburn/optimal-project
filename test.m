%% Test
clc
clear variables
close all

%% NCAT RUN DATA

data = load("4_10_ncat_lateral_3_2023-04-10-16-23-03.mat");

%% NOVATEL GPS

nov_lat = data.data.novatel_local.llhPositionTagged.latitude;
nov_lon = data.data.novatel_local.llhPositionTagged.longitude;
nov_alt = data.data.novatel_local.llhPositionTagged.altitude;

% lla from novatel
lla = [nov_lat; nov_lon; nov_alt];

% lla to ecef
ecef = lla2ecef(lla');

% lla to ned
xyzNED = lla2ned(lla', lla(:,1)', 'flat');

% ecef to ned
wgs84 = wgs84Ellipsoid('meter');

for i = 1:length(nov_lat)

[xN(i), yE(i), zd(i)] = ecef2ned(ecef(i,1), ecef(i,2), ecef(i,3), lla(1,1), lla(2,1), lla(3,1), wgs84);

end