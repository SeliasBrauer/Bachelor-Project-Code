clear variables; close all; clc; 

% Loading data.
load 'VORTALL'

% Creating data matrix . 
X = VORTALL;

%% Plot Vorticity
nx = 199;  % Number of grid points in x-direction
ny = 449;  % Number of grid points in y-direction
timestep = 69; % Can be chosen between 1 and 150. 
plotCylinder(reshape(X(:,timestep),nx,ny));
title('Vorticity Snapshot')

%% Compute POD modes
Xavg = mean(X,2); % Mean subtracted data matrix is found

%compute eigenmodes on mean subtracted data
X_e = X - Xavg*ones(1,size(X,2));
[U,S,V] = svd(X_e,'econ'); 

%plot first principal components
for i = 1:2
plotCylinder(reshape(U(:,i) ,nx,ny));
title(sprintf('Mode %i',i))
end
plotCylinder(reshape(Xavg ,nx,ny));
title('Mean')
%% Coefficient of modes time series
m = 2; %number of time series used
s = diag(S); %singular values vector
a = V.*s'; %coefficeint of modes time series. columns are times series for each mode 

a = a(:,1:m); %coefficients used further on

%figure()
%plot(a(:,1),a(:,2))


%% Finding derivatives of system amplitudes

dt = 0.2; % temp value, find real value
tspan = 0:dt:150*dt;

for i = 2:length(a)-1
    da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
end

figure() %plot coeficcients and derivatives
   subplot(2,1,1); hold on; grid on; 
   plot(1:149, da(:,1),"--b")
   plot(1:149, a(2:end-1,1),"-b")
   ylabel('Mode 1')

   subplot(2,1,2); hold on; grid on; 
   plot(1:149, da(:,2),"--r")
   plot(1:149, a(2:end-1,2),"-r")
   ylabel('Mode 2')
   
   xlabel('Time'); 

%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 3;
n = 2; %number of independent variables in system 
Theta = poolData(a,n,polyorder);
m = size(Theta,2);


%% Compute Sparse regression: sequential least squares

lambda = 0.1; % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta(2:end-1,:),da,lambda,n);
poolDataLIST({'x','y'},Xi,n,polyorder);

%% Compute antiderivative from sparse dynamics
x0 = a(1,:);
tspan = 0: dt: 150*dt;


%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
%[t,ai] = ode45(@(t,x) system(t,x,Xi), tspan, x0); 
[t,ai] = ode45(@(t,x) Diffeq_id_sys(t,x,Xi,n,polyorder), tspan, x0);

figure() %plot amplitudes of modes
subplot(2,1,1); hold on; 
plot(ai(:,1),'b')
plot(a(:,1),'r')
ylabel('Mode 1');
legend('Identified system','Full system')
title('System amplitude')

subplot(2,1,2); hold on
plot(ai(:,2),'b'); 
plot(a(:,2),'r')
xlabel('Time'); ylabel('Mode 2')
legend('Identified system','Full system')


figure() %plot phase space 
plot(ai(:,1), ai(:,2),'b');  grid on; hold on;
plot(a(:,1),a(:,2),'r',LineWidth=2); 
legend('Identified system','Full system')
xlabel('Mode 1'); ylabel('Mode 2')

%% Recreate flow from identified system and POD modes 

timestep = 69; % choose from 1 to 151

Recreate = U(:,1)*ai(timestep,1) + U(:,2)*ai(timestep,2) + Xavg; 

plotCylinder(reshape(Recreate,nx,ny)); 
title('Recreated vorticity Snapshot')

%% System of identified differential equations
function dadt = system(t,x,Xi)
g =[1; x(1); x(2); x(1)^2; x(1)*x(2); x(2)^2; x(1)^3; x(1)^2*x(2); x(1)*x(2)^2 ; x(2)^3];

dadt = (Xi.*g)';
dadt = [sum(dadt(1,:)); sum(dadt(2,:))];
end 
