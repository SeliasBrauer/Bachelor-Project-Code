clc; clear variables; close all; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


%load data unconfined: dt = 0.02;  5*dt = 0.1
dx = 1/22;
load VORTALL_unconfined_SINDy.mat
Lift_uncon = readmatrix("C:\Users\selia\Desktop\Unconfined mesh 22 presentation\Lift_unconfined_mesh_22.csv");

X = VORTALL_unconfined_SINDy;

Lift_uncon = Lift_uncon(:,2);

%cutoff data not included in vorticity data: 
idx = 5005 : 5 : length(Lift_uncon); 
Lift_uncon = Lift_uncon(idx); 

[~,I_max] = findpeaks(Lift_uncon); 

%cutoff vorticity data so that only whole periods are considered. 
Lift_uncon = Lift_uncon(I_max(1):I_max(end)-1);
X = X(:,I_max(1):I_max(end)-1);

vg = reshape(X,89,199,[]);
AnimateSquare(vg,'Unconfined_SINDy.gif',1/22,0.02,1);

%%
%load confined data dx = 1/30; dt = 0.01; 5*dt = 0.05
load VORTALL_confined_sindy.mat

X_con = VORTALL_confined_sindy;

Lift_con = readmatrix("C:\Users\selia\Desktop\Simulation data\Confined mesh 30\Lift_Mesh_30_confined.csv");

Lift_con = Lift_con(:,2);

%cutoff data not included in vorticity data: 
idx = 10005 : 5 : length(Lift_con); 
Lift_con = Lift_con(idx); 

%Cutoff half of vorticity data making 10*dt = 0.1
X_con = X_con(:,1:2:end);
Lift_con = Lift_con(1:2:end);

%find peak indexes. 
[~,I_max] = findpeaks(Lift_con); 

%cutoff vorticity data so that only whole periods are considered. 
Lift_con = Lift_con(I_max(1):I_max(end));
plot(Lift_con);
X_con = X_con(:,I_max(1):I_max(end)-1);

X_con = reshape(X_con, 121,271,[]);



AnimateSquare(X_con,'Confined_SINDy.gif',1/30,0.02,1);
