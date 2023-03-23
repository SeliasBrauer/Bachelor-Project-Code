clear variables; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Loading data.
load 'VORTALL'

% Creating data matrix . 
X = VORTALL;

%% Plot Vorticity
nx = 199;  % Number of grid points in x-direction
ny = 449;  % Number of grid points in y-direction
timestep = 10; % Can be chosen between 1 and 150. 
plotCylinder(reshape(X(:,timestep),nx,ny));
title('Vorticity Snapshot')

%% Compute POD modes
Xavg = mean(X,2); % Mean subtracted data matrix is found

%compute eigenmodes on mean subtracted data
X_B = X - Xavg*ones(1,size(X,2));
[U,S,V] = svd(X_B,'econ'); 

%plot first principal components
for i = 1:3
plotCylinder(reshape(U(:,i) ,nx,ny));
title(sprintf('Mode %i',i))
end
plotCylinder(reshape(Xavg ,nx,ny));
title('Mean')
%% Coefficient of modes time series
m = 3; %number of time series used
s = diag(S); %singular values vector

a = V.*s'; %coefficeint of modes time series. columns are times series for each mode 

a = a(:,1:m); %coefficients used further on

%% Finding derivatives of system amplitudes

%dt = 0.1984; % Value found using von karman

dt = 0.2; % real value of data.


tspan = 0:dt:150*dt;

 %compute derivative using finite difference
for i = 2:length(a)-1
    da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
end

figure() %plot coeficcients and derivatives
   subplot(3,1,1); hold on; 
   plot(1:149, da(:,1),"--b")
   plot(1:149, a(2:end-1,1),"-b")
   ylabel('Mode 1')
   
   subplot(3,1,2); hold on; 
   plot(1:149, da(:,2),"--r")
   plot(1:149, a(2:end-1,2),"-r")
   ylabel('Mode 2')

   subplot(3,1,3); hold on; grid on; 
   plot(1:149, da(:,3),"--g")
   plot(1:149, a(2:end-1,3),"-g")
   ylabel('Mode 3')

   xlabel('Time'); 


%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 3;
nVars = 3; %number of independent variables in system 
Theta = poolData(a,nVars,polyorder);
%m = size(Theta,2);


%% Compute Sparse regression: sequential least squares

lambda = 0.01; % lambda is our sparsification knob.

%find best fit coefficients
Xi = sparsifyDynamics(Theta(2:end-1,:),da,lambda,nVars);
poolDataLIST({'x','y','z'},Xi,nVars,polyorder);

%% Compute antiderivative from sparse dynamics
x0 = a(1,:); %initial values taken from time series amplitude

%extrapolate system in time to see if unstable
tspan = 0: dt: 300*dt;

%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
%integrate discovered dynamics with ode 45

[t,ai] = ode45(@(t,x) Diffeq_id_sys(t,x,Xi,nVars,polyorder), tspan, x0); 

%plot amplitudes of modes along with discovered amplitudes
figure() 
for i = 1:m
    subplot(m,1,i); hold on; 
    plot(ai(:,i),'b')
    plot(a(:,i),'r')
    ylabel(sprintf('Mode %i',i));
    legend('Identified system','Full system')

    if i == 1; title('System amplitude');
    elseif i == m; xlabel('Timestep');
    end
end

figure() %plot phase space with discovered
plot3(ai(:,1), ai(:,2),ai(:,3),'b');  grid on; hold on;
plot3(a(:,1),a(:,2),a(:,3),'r',LineWidth=2); 
legend('Identified system','Full system'); grid on ; 
xlabel('Mode 1'); ylabel('Mode 2'); zlabel('Mode 3')

%% Recreate flow from identified system and POD modes 

%timestep = 69; % choose from 1 to 151
%Sum all modes with discovered amplitudes. 
Recreate(:,:) = Xavg.*ones(size(VORTALL,1),length(tspan)); 
for j = 1:size(Recreate,2)
    for i = 1:m 
    Recreate(:,j) = Recreate(:,j) + U(:,i)*ai(j,i);
    end
end

%stack recreated data: 
for i = 1:size(Recreate,2)
    VORT_stack(:,:,i) = reshape(Recreate(:,i),nx,ny);
end

%AnimateCylinder(VORT_stack)
%title('Recreated vorticity Snapshot')

%% Energy diagram of first 10 modes: 

%plot the energy camptured by first 10 modes
for i = 1:10
energy_cap(i) = sum(s(1:i)) / sum(s) * 100; 
end 

figure(); 
plot(energy_cap,'or',MarkerFaceColor='red'); grid on; 
title('Energy caputed by the first 10 modes')
xlabel('Number of modes included')
ylabel('Energy captured by the model (\%)')
xlim([0 10]); ylim([0 100]); 

ME_CYL = energy_cap; 
save("Mode_energy","ME_CYL",'-append');
