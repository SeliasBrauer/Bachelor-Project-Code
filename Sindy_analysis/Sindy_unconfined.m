clear variables; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

m = 6; %number of modes used

%force symmetry? yes = 1; no = 0; 
symmetry = 1;

% Loading data.
load VORTALL_UNCONFINED_LY10.mat

stepsize = 1; %define data sampling rate fx: stepsize = 2 : every second snapshot is used.
% Creating data matrix . 
X = VORTALL_UNCONFINED_LY10(:,1:stepsize:end);

clear VORTALL_UNCONFINED_LY10

%grid size
dt = 0.02*stepsize; % real value of data. dt = 0.01;
dx = 1/22; 
nx = 10  *22 + 1;  % Number of grid points in x-direction
ny = 199;  % Number of grid points in y-direction

%% Plot Vorticity

timestep = 10; % Can be chosen between 1 and 217. 
plotSquare(reshape(X(:,timestep),nx,ny),dx);
title('Vorticity Snapshot')

%% Compute POD modes
if symmetry == 1 %Creating symmetrized data matrix for POD.

    Y = [X X]; % Y will be the symmetrized data matrix.
    for k=1:size(X,2)

        % Flipping y-coordinate.
        xflip = reshape(flipud(reshape(X(:,k),nx,ny)),ny*nx,1);

        % Adding sign change.
        Y(:,k+size(X,2)) = -xflip;
    end

    X = Y;
    plotSquare(reshape(X(:,1),nx,ny),dx); % plot of wake.

    plotSquare(reshape(X(:,1+size(X,2)/2),nx,ny),dx); % plot of transformed wake.
   
end

Xavg = mean(X,2); % Mean subtracted data matrix is found

%compute eigenmodes on mean subtracted data
X_B = X - Xavg*ones(1,size(X,2));
[U,S,V] = svd(X_B,'econ'); 

%plot first principal components
for i = 1:m
subplot(m/2, 2,i)
plotSquare(reshape(U(:,i) ,nx,ny),dx);
title(sprintf('Mode %i',i))
end

figure()
plotSquare(reshape(Xavg ,nx,ny),dx);
title('Mean')
%% Coefficient of modes time series
s = diag(S); %singular values vector

a = V.*s'; %coefficeint of modes time series. columns are times series for each mode 

a = a(:,1:m); %coefficients used further on

for i = 1:m 
    [~,pktimes] = findpeaks(a(:,i));
    Period(i) = mean(diff(pktimes))*dt;
end
freq = 1./Period;

plot(freq,'ok','MarkerFaceColor','k'); grid on; 
xlabel('Mode'); ylabel('Freqeuncy [Hz]')

%% Finding derivatives of system amplitudes

if symmetry == 1 
    tspan = 0:dt:(size(X,2)/2 - 1)*dt;

    for i = 2:size(X,2)/2 - 1
        da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    for i = size(X,2)/2 + 2  : size(X,2) - 1
        da(i-3,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    figure() %plot coeficcients and derivatives 
    for i = 1:m
        subplot(m,1,i); hold on; 
        plot(tspan(2:end-1), da(1:end/2,i),".b")
        plot(tspan(2:end-1), a(2:end/2-1,i),".r")
        legend('Derivative','Amplitude')
        ylabel(sprintf('Mode %i',i))
    end
    xlabel('Time'); 

else
    tspan = 0:dt:(size(X,2) - 1)*dt;

    %compute derivative using finite difference
    for i = 2:length(a)-1
        da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    figure() %plot coeficcients and derivatives 
    for i = 1:m
        subplot(m,1,i); hold on; 
        plot(tspan(2:end-1), da(:,i),"--b")
        plot(tspan(2:end-1), a(2:end-1,i),"-r")
        legend('Derivative','Amplitude')
        ylabel(sprintf('Mode %i',i))
    end
    xlabel('Time'); 
end

% save derivative and time span for convergance test of finite difference. 
struct_020.dt = dt; 
struct_020.timeseries = tspan;
struct_020.derivative = da;
%save("Derivatives_unconfined", "struct_020",'-append');


%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 2;
nVars = m; %number of independent variables in system 
Theta = poolData_nconstant(a,nVars,polyorder);

%% Compute Sparse regression: sequential least squares

lambda = 0.00031; 
lambda2 = [0.5 , 0.5, 0.0075, 0.0075, 0.04, 0.08]; % lambda is our sparsification knob.
lambda = 0.00031; 
%lambda2 = [0,0,0,0,0,0]; 
C = []; d = []; 
%constrainst for 4 modes
%{
C{1} = [3,4,1];  d = 0;
C{2} = [3,8,1];  d = [d;0];
C{3} = [3,11,1]; d = [d;0];
C{4} = [3,13,1]; d = [d;0];
C{5} = [3,14,1]; d = [d;0];

C{6} = [4,3,1];  d = [d;0];
C{7} = [4,7,1];  d = [d;0];
C{8} = [4,10,1]; d = [d;0];
C{9} = [4,12,1]; d = [d;0];
C{10} = [4,13,1]; d = [d;0];
%}


%constraints for 6 modes
%z not connected to alpha, beta or gamma
C{1} = [3,4,1];  d = 0;
C{2} = [3,5,1];  d = [d;0];
C{3} = [3,6,1];  d = [d;0];
C{4} = [3,10,1]; d = [d;0];
C{5} = [3,11,1]; d = [d;0];
C{6} = [3,12,1]; d = [d;0];
C{7} = [3,15,1]; d = [d;0];
C{8} = [3,16,1]; d = [d;0];
C{9} = [3,17,1]; d = [d;0];
C{10} = [3,19,1]; d = [d;0];
C{11} = [3,20,1]; d = [d;0];
C{12} = [3,21,1]; d = [d;0];
C{13} = [3,22,1]; d = [d;0];
C{14} = [3,23,1]; d = [d;0];
C{15} = [3,24,1]; d = [d;0];
C{16} = [3,25,1]; d = [d;0];
C{17} = [3,26,1]; d = [d;0];
C{18} = [3,27,1]; d = [d;0];

%alpha not connected to z, beta or gamma
C{19} = [4,3,1];  d = [d;0];
C{20} = [4,5,1];  d = [d;0];
C{21} = [4,6,1];  d = [d;0];
C{22} = [4,9,1];  d = [d;0];
C{23} = [4,11,1];  d = [d;0];
C{24} = [4,12,1];  d = [d;0];
C{25} = [4,14,1]; d = [d;0];
C{26} = [4,16,1]; d = [d;0];
C{27} = [4,17,1]; d = [d;0];
C{28} = [4,18,1]; d = [d;0];
C{29} = [4,19,1]; d = [d;0];
C{30} = [4,20,1]; d = [d;0];
C{31} = [4,21,1]; d = [d;0];
C{32} = [4,23,1]; d = [d;0];
C{33} = [4,24,1]; d = [d;0];
C{34} = [4,25,1]; d = [d;0];
C{35} = [4,26,1]; d = [d;0];
C{36} = [4,27,1]; d = [d;0];

%beta not connected to gamma
C{37} = [5,6,1];  d = [d;0];
C{38} = [5,12,1]; d = [d;0];
C{39} = [5,17,1]; d = [d;0];
C{40} = [5,21,1]; d = [d;0];
C{41} = [5,24,1]; d = [d;0];
C{42} = [5,26,1]; d = [d;0];
C{43} = [5,27,1]; d = [d;0];

%gamma not connected to beta
C{44} = [6,5,1];  d = [d;0];
C{45} = [6,11,1]; d = [d;0];
C{46} = [6,16,1]; d = [d;0];
C{47} = [6,20,1]; d = [d;0];
C{48} = [6,23,1]; d = [d;0];
C{49} = [6,25,1]; d = [d;0];
C{50} = [6,26,1]; d = [d;0];

%x not connected to z, alpha, beta or gamma
C{51} = [1,3,1]; d = [d;0];
C{52} = [1,4,1]; d = [d;0];
C{53} = [1,5,1]; d = [d;0];
C{54} = [1,6,1]; d = [d;0];
C{55} = [1,9,1]; d = [d;0];
C{56} = [1,10,1]; d = [d;0];
C{57} = [1,11,1]; d = [d;0];
C{58} = [1,12,1]; d = [d;0];
C{59} = [1,14,1]; d = [d;0];
C{60} = [1,15,1]; d = [d;0];
C{61} = [1,16,1]; d = [d;0];
C{62} = [1,17,1]; d = [d;0];
C{63} = [1,18,1]; d = [d;0];
C{64} = [1,19,1]; d = [d;0];
C{65} = [1,20,1]; d = [d;0];
C{66} = [1,21,1]; d = [d;0];
C{67} = [1,22,1]; d = [d;0];
C{68} = [1,23,1]; d = [d;0];
C{69} = [1,24,1]; d = [d;0];
C{70} = [1,25,1]; d = [d;0];
C{71} = [1,26,1]; d = [d;0];
C{72} = [1,27,1]; d = [d;0];

%y not connected to z, alpha, beta or gamma
C{73} = [2,3,1]; d = [d;0];
C{74} = [2,4,1]; d = [d;0];
C{75} = [2,5,1]; d = [d;0];
C{76} = [2,6,1]; d = [d;0];
C{77} = [2,9,1]; d = [d;0];
C{78} = [2,10,1]; d = [d;0];
C{79} = [2,11,1]; d = [d;0];
C{80} = [2,12,1]; d = [d;0];
C{81} = [2,14,1]; d = [d;0];
C{82} = [2,15,1]; d = [d;0];
C{83} = [2,16,1]; d = [d;0];
C{84} = [2,17,1]; d = [d;0];
C{85} = [2,18,1]; d = [d;0];
C{86} = [2,19,1]; d = [d;0];
C{87} = [2,20,1]; d = [d;0];
C{88} = [2,21,1]; d = [d;0];
C{89} = [2,22,1]; d = [d;0];
C{90} = [2,23,1]; d = [d;0];
C{91} = [2,24,1]; d = [d;0];
C{92} = [2,25,1]; d = [d;0];
C{93} = [2,26,1]; d = [d;0];
C{94} = [2,27,1]; d = [d;0];
%}

%find best fit coefficients
if symmetry == 1
    %definde used values corrosponding to 
    range = [2:(size(X,2)/2-1) , (size(X,2)/2 +2) : (size(X,2)-1)];
    
    %use constrained sparsifying function

    %Xi = sparsifyDynamics(Theta(range,:),da,lambda,nVars);
    
    %Xi = sparsifyDynamics_con(Theta(range,:),da,lambda2,nVars,C,d);

    Xi = sparsifyDynamics_con_single(Theta(range,:),da, lambda2, nVars,C,d);

    %Xi = sparsifyDynamics_con_lasso(Theta(range,:),da,lambda,nVars,C,d);

    %Xi = sparsifyDynamics_con_mix(Theta(range,:),da,lambda,lambda2,nVars,C,d);
else
    Xi = sparsifyDynamics(Theta(2:end-1,:),da,lambda,nVars);
end


% list of variable names
var_name = {'x','y','z','alpha','beta','gamma','i','j','k','v'};
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
if symmetry ==1 
    Plt_inx = 1:size(a,1)/2;
else
    Plt_inx = 1:size(a,1);
end

for i = 1:m
    subplot(m,1,i); hold on; grid on; 
    plot(t,ai(:,i),'b')
    plot(t(Plt_inx) , a(Plt_inx,i),'r')
    ylabel(sprintf('Mode %i',i));
    legend('Identified system','Full system')

    if i == 1; title('System amplitude');
    elseif i == m; xlabel('Time [s]');
    end
end


%plot phase space with discovered dynamics
figure()
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
        plot(ai(:,i), ai(:,j),'b');  grid on; hold on; axis equal 
        plot(a(Plt_inx,i),a(Plt_inx,j),'r',LineWidth=2); 
        %legend('Identified system','Full system'); grid on ; 
        xlabel(sprintf('Mode %i',i)); ylabel(sprintf('Mode %i',j));
        end
    end 
end 


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


ME_SQ_UNCON = energy_cap; 
%save("Mode_energy","ME_SQ_UNCON",'-append');
