clc; clear variables; close all; 

%import vorticity data for confined flow around square cylinder with Ly = 5;
load VORTALL_CONFINED_LY5.mat

X = VORTALL_CONFINED_LY5;

% grid properties. 
dx = 1/22;
nx = 89; ny = 199; 

%import lift coefficient corrospondign to vorticity data
Lift = readmatrix('C:\Users\selia\Desktop\Vorticity confined, Ly = 5\Lift_Confined_LY5.csv');
Lift = Lift(:,2);

%cutoff data not included in vorticity data: 
idx = 5010 : 10 : length(Lift); 
Lift = Lift(idx); 

[~,I_max] = findpeaks(Lift);

%cutoff vorticity data so that only whole periods are considered. 
Lift = Lift(I_max(1):I_max(end));
VORTALL_CONFINED_LY5_SYMM = X(:,I_max(1):I_max(end));

%save modified data for SINDy analysis. 
save("VORTALL_CONFINED_LY5_SYMM","VORTALL_CONFINED_LY5_SYMM");

%create stacked matix for plotting: 
X_stack = reshape(VORTALL_CONFINED_LY5_SYMM,nx,ny,[]);

%AnimateSquare(X_stack,'Square_confined_LY5_symmetric_anim.gif',0.05,1);
