clear variables; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

m = 6; %number of modes used

%force symmetry? yes = 1; no = 0; 
symmetry = 1;

% Loading data.
load VORTALL_SQUARE_UNCONFINED.mat

% Creating data matrix . 
X = VORTALL_SQUARE_UNCONFINED;

%grid size
dx = 1/22; 

%% Plot Vorticity
nx = 89;  % Number of grid points in x-direction
ny = 199;  % Number of grid points in y-direction

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
for i = 1:2
plotSquare(reshape(U(:,i) ,nx,ny),dx);
title(sprintf('Mode %i',i))
end
   

plotSquare(reshape(Xavg ,nx,ny),dx);
title('Mean')
%% Coefficient of modes time series
s = diag(S); %singular values vector

%Dirty fix making mode energies equal. 
%{
s(1:2) = mean(s(1:2)); 
s(3:4) = mean(s(3:4));
s(5:6) = mean(s(5:6));
%}


a = V.*s'; %coefficeint of modes time series. columns are times series for each mode 

a = a(:,1:m); %coefficients used further on

%% Finding derivatives of system amplitudes

dt = 0.2; % real value of data.

if symmetry == 1 
    tspan = 0:dt:(size(X,2)/2 - 1)*dt;

    for i = 2:size(X,2)/2 - 1
        da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    for i = size(X,2)/2 + 2  : size(X,2) - 1
        da(i-3,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

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

%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 2;
nVars = m; %number of independent variables in system 
Theta = poolData_nconstant(a,nVars,polyorder);

%% Compute Sparse regression: sequential least squares

lambda = 0.2; % lambda is our sparsification knob.

%find best fit coefficients
if symmetry == 1
    range = [2:(size(X,2)/2-1) , (size(X,2)/2 +2) : (size(X,2)-1)];

    %Xi = sparsifyDynamics(Theta(range,:),da,lambda,nVars);
    Xi = sparsifyDynamics_con(Theta(range,:),da,lambda,nVars,[],[]);

    %Xi = lsqlin(Theta(range,:), da, [],[]);
    % use lasso regression instead of STLS
    %{
    for i = 1:m
    [Xi(:,i),stats] = lasso(Theta(range,:),da(:,i),'Lambda',lambda);
    end
    %}
else
    Xi = sparsifyDynamics(Theta(2:end-1,:),da,lambda,nVars);
end

%%force constraints
%{
m1 = 1; %( abs(Xi(2,1)) + abs(Xi(1,2)) ) / 2;
Xi(2,1) = - m1 ;  
Xi(1,2) = m1; 

m2 = 2; % ( abs(Xi(4,3)) + abs(Xi(3,4)) ) / 2;
Xi(4,3) =  m2 ;  
Xi(3,4) = - m2 ;  

m3 = 3; %( abs(Xi(6,5)) + abs(Xi(5,6)) ) / 2;
Xi(6,5) = - m3 ;  
Xi(5,6) =  m3 ;  
%}


% list of variable names
var_name = {'x','y','z','alpha','beta','gamma','i','j','k','v'};
%print dynamics
poolDataLIST_nconstant(var_name(1:m),Xi,nVars,polyorder);

%% Compute antiderivative from sparse dynamics
x0 = a(1,:); %initial values taken from time series amplitude

%extrapolate system in time to see if unstable
tspan = 0: dt: 217*dt + 10*31*dt;

%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
%integrate discovered dynamics with ode 45

[t,ai] = ode45(@(t,x) Diffeq_id_sys_nconstant(t,x,Xi,nVars,polyorder), tspan, x0); 

%plot amplitudes of modes along with discovered amplitudes
figure() 
for i = 1:m
    subplot(m,1,i); hold on; 
    plot(t,ai(:,i),'b')
    plot(t(1:length(a(1:end/2,i))) , a(1:end/2,i),'r')
    ylabel(sprintf('Mode %i',i));
    legend('Identified system','Full system')

    if i == 1; title('System amplitude');
    elseif i == m; xlabel('Time [s]');
    end
end

%plot phase space with discovered dynamics
figure()

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
        plot(a(1:end/2,i),a(1:end/2,j),'r',LineWidth=2); 
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
