clc; clear variables; close all; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%load data unconfined: dt = 0.02;  5*dt = 0.1
load VORTALL_unconfined_presentation.mat
Lift_con = readmatrix("C:\Users\selia\Desktop\Unconfined mesh 22 presentation\Unconfined_22_lift_force.csv");

X = VORTALL_unconfined_presentation;

Lift_con = Lift_con(:,2);

%cutoff data not included in vorticity data: 
idx = 5005 : 5 : length(Lift_con); 
Lift_con = Lift_con(idx); 

[~,I_max] = findpeaks(Lift_con); 

%cutoff vorticity data so that only whole periods are considered. 
Lift_con = Lift_con(I_max(1):I_max(end));
X = X(:,I_max(1):I_max(end)-1);

vg = reshape(X,89,683,[]);
AnimateSquare(vg,'Unconfined_presentation.gif',1/22,0.02,1);

%%
%load confined data dx = 1/30; dt = 0.01; 5*dt = 0.05
load VORTALL_confined_presentation.mat
X_con = VORTALL_confined_presentation;

Lift_con = readmatrix("C:\Users\selia\Desktop\Unconfined mesh 22 presentation\Confined_30_lift_force.csv");

Lift_con = Lift_con(:,2);

%cutoff data not included in vorticity data: 
idx = 10005 : 10 : length(Lift_con); 
Lift_con = Lift_con(idx); 

%Cutoff half of vorticity data making 10*dt = 0.1
X_con = X_con(:,:,1:2:end);

%find peak indexes. 
[~,I_max] = findpeaks(Lift_con); 

%cutoff vorticity data so that only whole periods are considered. 
Lift_con = Lift_con(I_max(1):I_max(end));
X_con = X_con(:,I_max(1):I_max(end)-1);

X_con = reshape(X, 121,931,[]);



AnimateSquare(X_con,'Confined_presentation.gif',1/30,0.02,1);
