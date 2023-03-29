clc; close all; clear variables
%% Finding the timestep size from vorticity data: 

% loading vorticity data
load 'VORTALL'
load 'VELOCITIES.mat'
X = VORTALL; 
nx = 199; % # grid points in x
ny = 449; % # grid points in y

%plot snapshot of data
timestep = 1; % choose between 1 and 151

%plotCylinder(reshape(X(:,timestep),nx,ny))



%% Find number of timesteps for a period
for i = 1: 151
   %plot error between chosen timestep and the rest
error(i) = sum(X(:,timestep) - X(:,i)) ;
end 

findpeaks(error) % find peaks in error
title('Von Karman frequency')
ylabel('Error'); xlabel('Timestep')
[pks,locs] = findpeaks(error);

% find number of timestep for a period: 
timestep_period = mean(locs(2:end) - locs(1:end-1)); 

%% Find bulk velocity along the stream for data
plotCylinder(reshape(UALL(:,timestep),nx,ny))
plotCylinder(reshape(VALL(:,timestep),nx,ny))

Uavg = mean(mean(UALL));  %Bulk velocity mean is found: 
%% Find theoretical Von Karman frequency for Re = 100
% standard porpperties of air and water are found to determine the
% most likely diamenter for the cylinder.
% Re = rho*D*Uavg / mu ; 
% frequency found as  f = 0.21 * V/D (1 - 20/Re) 

U = 1; 
D = 1; 
nu = 0.01;

rho_w = 998; % kg/m^3 
mu_w  = 1.002*10^-3; %Pa*s

rho_a = 1.204; %kg/m^3 
mu_a = 1.825*10^-5; %Pa*s

%diameters are determined for the cylinder. D = Re*mu / (rho*Uavg) 

D_w = 100*mu_w / (rho_w*Uavg);
D_a = 100*mu_a / (rho_a*Uavg);



VKF = 0.21* U/D*(1-20/100);
VKF_w = 0.21* Uavg/D_w*(1-20/100); % 0.21 is for range  Re >= 1000
VKF_a = 0.21* Uavg/D_a*(1-20/100); % 0.21 is for range  Re >= 1000

T = 1/VKF;
T_w = 1/VKF_w; %time period
T_a = 1/VKF_a; %time period

dt = T/timestep_period
dt_w = T_w/timestep_period; 
dt_a = T_a/timestep_period; 



