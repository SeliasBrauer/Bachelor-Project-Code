clear variables; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

m = 6; %number of modes used

% Loading data.
load VELOCITIES_circ.mat

% Creating data matrix . 
X = [UALL;VALL]; % stack velocity data

%index splitting data: 
n_half = height(X)/2;

nx = 199;  % Number of grid points in x-direction
ny = 449;  % Number of grid points in y-direction

%% Plot Vorticity
timestep = 10; % Can be chosen between 1 and 217. 

plotCylinder(reshape(X(1:n_half,timestep),nx,ny));
title('Velocity-x Snapshot')

plotCylinder(reshape(X(n_half+1:end,timestep),nx,ny));
title('Velocity-y Snapshot')

%% Compute POD modes 

Xavg = mean(X,2); % Mean subtracted data matrix is found

%compute eigenmodes on mean subtracted data
X_B = X - Xavg*ones(1,size(X,2));

[U,S,V] = svd(X_B,'econ'); 

%flip signs such that modes are consistent. 
for i = 1:size(S,1)
    U(:,i) = U(:,i) * V(1,i)/abs(V(1,i));
    V(:,i) = V(:,i) * V(1,i)/abs(V(1,i));
end


%% plot first POD modes normalized. 
Xavg_norm = Xavg(1:n_half) ./ max(max(abs( Xavg(1:n_half) )));
Yavg_norm = Xavg(n_half+1:end) ./ max(max(abs(Xavg(n_half+1:end)))); 

for i = 1:m
    Ux_norm(1:n_half,i) = U(1:n_half,i) ./ max(max(abs(U(1:n_half,i))));
    Uy_norm(:,i) = U(n_half+1:end,i) ./ max(max(abs(U(n_half+1:end,i))));
end


mode_fig = figure('Position',[400,10,300,700]); 
subplot(m+1,1,1)

plotCylinder(reshape(Xavg_norm ,nx,ny));
title('Mean')

for i = 1:m
subplot(m+1, 1,i+1)
plotCylinder(reshape(Ux_norm(:,i) ,nx,ny));
title(sprintf('Mode %i',i))

end


savepath = "C:\Users\selia\Desktop\Bachelor Project images\Qualitative analysis"; 
mode_fig_name = "Velocity-x modes.pdf"; 
exportgraphics(mode_fig, fullfile(savepath,mode_fig_name), "ContentType","vector");

mode_fig = figure('Position',[400,10,300,700]); 
subplot(m+1,1,1)

plotCylinder(reshape(Yavg_norm ,nx,ny));
title('Mean')

for i = 1:m
subplot(m+1, 1,i+1)
plotCylinder(reshape(Uy_norm(:,i) ,nx,ny));
title(sprintf('Mode %i',i))
end


savepath = "C:\Users\selia\Desktop\Bachelor Project images\Qualitative analysis"; 
mode_fig_name = "Velocity-y modes.pdf"; 
exportgraphics(mode_fig, fullfile(savepath,mode_fig_name), "ContentType","vector");
%% Coefficient of modes time series
s = diag(S); %singular values vector

a = V.*s'; %coefficeint of modes time series. columns are times series for each mode 

a = a(:,1:m); %coefficients used further on

%% Finding derivatives of system amplitudes

dt = 0.2; % real value of data.
tspan = 0:dt:(size(X,2) - 1)*dt;

%compute derivative using finite difference
for i = 2:length(a)-1
    da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
end
    
amp_fig = figure('position',[400,100,800,500]); grid on; %plot coeficcients and derivatives 
for i = 1:m
    subplot(m,1,i); hold on; 
    plot(tspan(2:end-1), da(:,i),".b")
    plot(tspan(2:end-1), a(2:end-1,i),".r")
    if i == 1; legend('Derivative','Amplitude'); end
    ylabel(sprintf('$ a_{%i}$',i))
    set(gca,'fontsize',10); grid on; 
end
set(gcf,'color','w');
xlabel('Time'); 

savepath = "C:\Users\selia\Desktop\Bachelor Project images\Qualitative analysis"; 
amp_fig_name = "Velocity amplitudes.pdf"; 
%exportgraphics(amp_fig, fullfile(savepath,amp_fig_name), "ContentType","vector");
%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 2;

nVars = m; %number of independent variables in system 

%pool data adding velocity modes
Theta = poolData_nconstant(a,nVars,polyorder);

%% Compute Sparse regression: sequential least squares

lambda = 0.2; % lambda is our sparsification knob.


Xi = sparsifyDynamics(Theta(2:end-1,:),da,lambda,nVars);


% list of variable names
var_name = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','a11','a12','a13','a14','a15','a16','a17','a18','a19','a20'};
%print dynamics
poolDataLIST_nconstant(var_name(1:m),Xi,nVars,polyorder);

%% Compute antiderivative from sparse dynamics

x0 = a(1,:); %initial values taken from time series amplitude

%extrapolate system in time to see if unstable
tspan = 0: dt: 1 * size(a,1)*dt - dt;

%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

%integrate discovered dynamics with ode 45

[t,ai] = ode45(@(t,x) Diffeq_id_sys_nconstant(t,x,Xi,nVars,polyorder), tspan, x0); 

%plot amplitudes of modes along with discovered amplitudes
figure()
set(gcf,'color','white');
%define plotting indexes
Plt_inx = 1:size(a,1);


for i = 1:m
    subplot(m,1,i); hold on; grid on; 
    plot(t,ai(:,i),'b')
    plot(t(Plt_inx) , a(Plt_inx,i),'r')
    ylabel(sprintf('$a_{%i}$',i));
    
    if i == 1; legend('Reproduced','Original');
    elseif i == m; xlabel('Time [s]');
    end
end


%plot phase space with discovered dynamics
fig_phase = figure('position',[400,100,800,500]);
set(gcf,'color','white');
%vector of plot position very stupid way of finding vector of plot
%positions. 

pos = 1:(m-1)*(m-1);
pos = reshape(pos,m-1,m-1)';
for i = 1:m-1
    for j = 1:m-1
        if j> i 
            pos(j,i) = 0; 
        end 
    end
end
pos = nonzeros (reshape(pos',1,(m-1)*(m-1)));
k = 0; %counter
for i = 1:m 
    for j = 1:m
        if j > i
        k = k+1;
        subplot(m-1,m-1,pos(k))
        grid on; hold on; axis equal 
        plot(a(Plt_inx,i),a(Plt_inx,j),'r',LineWidth=2); 
        %plot(ai(:,i), ai(:,j),'b');
        if pos(k) == 1
        %lgd = legend('Full system','Identified system'); 
        %lgd.Position = [0.18,0.73,0.01,0.01]; 
        end
        xlabel(sprintf('$a_{%i}$',i)); ylabel(sprintf('$a_{%i}$',j));
        end
        set(gca,'fontsize',9)
    end 
end  

savepath = "C:\Users\selia\Desktop\Bachelor Project images\Qualitative analysis"; 
phase_fig_name = "Velocity Phase portrait.pdf"; 
exportgraphics(fig_phase, fullfile(savepath,phase_fig_name), "ContentType","vector");

%% Recreate flow from identified system and POD modes 


%Sum all modes with discovered amplitudes. 
Recreate(:,:) = Xavg.*ones(size(X,1),length(tspan)); 
for j = 1:size(Recreate,2)
    for i = 1:m 
    Recreate(:,j) = Recreate(:,j) + U(:,i)*ai(j,i);
    end
end

%stack recreated data: 
for i = 1:size(Recreate,2)
    VORT_stack(:,:,i) = reshape(Recreate(:,i),nx,ny);
end

%AnimateSquare(VORT_stack,'Square_flow_recreated_anim.gif',0.05,1);
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


ME_circ_vel = energy_cap; 
%uncomment to save energy data (default : leave as comment) 
%save("Mode_energy","ME_circ_vel",'-append');
