clc; clear variables; close all; 

%import vorticity data for unconfined flow around square cylinder;
load 'VORSTACK'
VORSTACK = vg; 

% grid properties. 
dx = 1/22;
nx = size(VORSTACK,1); ny = size(VORSTACK,2); nt = size(VORSTACK,3); 

%create Vortall matrix of vorticity data; 
for i = 1:nt
    X(:,i) = reshape(VORSTACK(:,:,i),1,nx*ny); 
end 

VORTALL_SQUARE_UNCONFINED = X;
save("VORTALL_SQUARE_UNCONFINED","VORTALL_SQUARE_UNCONFINED")

%import lift coefficient corrospondign to vorticity data
Lift = readmatrix('C:\Users\selia\Desktop\vorticity unconfined-20230303T215318Z-001\Lift_22x22.csv');
Lift = Lift(:,2);

%cutoff data not included in vorticity data: 
idx = 5010 : 10 : length(Lift); 
Lift = Lift(idx); 


[~,I_max] = findpeaks(Lift); 

%cutoff vorticity data so that only whole periods are considered. 
Lift = Lift(I_max(1):I_max(end));
VORTALL_SQUARE_UNCONFINED_SYMM = X(:,I_max(1):I_max(end));

%save modified data for SINDy analysis. 
save("VORTALL_SQUARE_UNCONFINED_SYMM","VORTALL_SQUARE_UNCONFINED_SYMM");

%create stacked matix for plotting: 
X_stack = reshape(VORTALL_SQUARE_UNCONFINED,nx,ny,[]);

%AnimateSquare(X_stack,0.05,1);
