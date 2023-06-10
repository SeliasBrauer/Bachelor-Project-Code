clear variables; close all; clc; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

m = 6; %number of modes used

%force symmetry? yes = 1; no = 0; 
symmetry = 1;

% Loading data.
load VORTALL_circ.mat

stepsize = 1; %define data sampling rate fx: stepsize = 2 : every second snapshot is used.
% Creating data matrix . 
X = - VORTALL_circ(:,1:stepsize:end);

clear VORTALL_circ;

%grid size
dt = 0.2*stepsize; % real value of data. dt = 0.01;
%dx = 1/22; 
ny = 449;  % Number of grid points in x-direction
nx = 199;  % Number of grid points in y-direction

%% Plot Vorticity
figure('Position',[400,500,2000,190])
timestep = 10; % Can be chosen between 1 and 217. 
plotCylinder(reshape(X(:,timestep),nx,ny));
set(gca,'fontsize',14)
%title('Vorticity Snapshot')

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
    
    figure('Position',[400,10,1000,200])
    subplot(1,2,1);
    plotCylinder(reshape(X(:,1),nx,ny)); % plot of wake.
    set(gca,'fontsize',14); title('Original'); 
    
    subplot(1,2,2); 
    plotCylinder(reshape(X(:,1+size(X,2)/2),nx,ny)); % plot of transformed wake.
    set(gca,'fontsize',14); title('Flipped');
end

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
Xavg_norm = Xavg ./ max(max(abs(Xavg)));

for i = 1:m
    U_norm(:,i) = U(:,i) ./ max(max(abs(U(:,i))));
end

mode_fig = figure('Position',[400,10,300,700]); 
subplot(m+1,1,1)

plotCylinder(reshape(Xavg_norm ,nx,ny));
title('Mean')

for i = 1:m
subplot(m+1, 1,i+1)
plotCylinder(reshape(U_norm(:,i) ,nx,ny));
title(sprintf('Mode %i',i))
end

savepath = "C:\Users\selia\Desktop\Bachelor Project images\Qualitative analysis"; 
mode_fig_name = "Circle modes.pdf"; 
%exportgraphics(mode_fig, fullfile(savepath,mode_fig_name), "ContentType","vector");

%% Coefficient of modes time series
s = diag(S); %singular values vector

a = V.*s'; %coefficeint of modes time series. columns are times series for each mode 

a = a(:,1:m); %coefficients used further on

for i = 1:m 
    if symmetry == 1
        [~,pktimes] = findpeaks(a(1:end/2,i));
        Period(i) = mean(diff(pktimes))*dt;
    else
        [~,pktimes] = findpeaks(a(:,i));
        Period(i) = mean(diff(pktimes))*dt;
    end
end
freq = 1./Period;

plot(freq,'ok','MarkerFaceColor','k'); grid on; 
xlabel('Mode'); ylabel('Freqeuncy [Hz]')

%% Finding derivatives of system amplitudes
div_fig = figure('position',[400,100,800,500]); grid on;
if symmetry == 1 
    tspan = 0:dt:(size(X,2)/2 - 1)*dt;

    for i = 2:size(X,2)/2 - 1
        da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    for i = size(X,2)/2 + 2  : size(X,2) - 1
        da(i-3,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    for i = 1:m
        subplot(m,1,i); hold on; grid on;
        plot(tspan(2:end-1), da(1:end/2,i),".b")
        plot(tspan(2:end-1), a(2:end/2-1,i),".r")
         if i == 1 
             legend('Derivative','Amplitude'); 
         end
        ylabel(sprintf('$a_{%i}$',i))
        set(gca,'fontsize',10)
    end
    xlabel('Time [s]'); 

else
    tspan = 0:dt:(size(X,2) - 1)*dt;

    %compute derivative using finite difference
    for i = 2:length(a)-1
        da(i-1,:) = (a(i+1,:) - a(i-1,:))./ (2*dt) ;
    end

    for i = 1:m
        subplot(m,1,i); hold on; grid on;
        plot(tspan(2:end-1), da(:,i),".b")
        plot(tspan(2:end-1), a(2:end-1,i),".r")
        legend('Derivative','Amplitude')
        ylabel(sprintf('Mode %i',i))
    end
    xlabel('Time [s]'); 
end
set(gcf,'color','w'); %background color white

% save derivative and time span for convergance test of finite difference. 
struct_020.dt = dt; 
struct_020.timeseries = tspan;
struct_020.derivative = da;
%save("Derivatives_unconfined", "struct_020",'-append');


savepath = "C:\Users\selia\Desktop\Bachelor Project images\Qualitative analysis"; 
div_fig_name = "Circle Amplitudes.pdf"; 
exportgraphics(div_fig, fullfile(savepath,div_fig_name), "ContentType","vector");


%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 2;
nVars = m; %number of independent variables in system 
Theta = poolData_nconstant(a,nVars,polyorder);

%% Compute Sparse regression: sequential least squares

lambda = 0.00031; 
lambda2 = [0.2 , 0.5, 0.0075, 0.004, 0.01, 0.006, 0.04, 0.02, 0.01, 0.0002]; % works well for 8 modes STLS

lambda2 = [0.2 , 0.5, 0.0075, 0.004, 0.02, 0.01, 0.01, 0.02, 0.01, 0.0002]; % works very well or 8 mode elastic net
lambda = 0.002; 
%alpha = 0.4; 


%lambda2 = [0.2 , 0.5, 0.005, 0.003, 0.002, 0.001, 0.01, 0.02, 0.01, 0.0002]; % works very well or 8 mode elastic net
%lambda = 0.08; 
%alpha = 0.3; 

%lambda2 = [0.2 , 0.5, 0.005, 0.003, 0.004, 0.001]; % works very well or 8 mode elastic net
lambda2 = 0.2* ones(1,m); 
%lambda = 0.02;
%alpha = 0.4; 



%lambda2 = [0.2 , 0.2, 0.001, 0.002, 0.005, 0.001]; % works very well or 8 mode elastic net
%lambda = 0.02;
%alpha = 1; 

C = []; d = [];

%[C,d] = hierarchical_con(nVars,2);
%C{end+1} = [4,8,1]; d = [d;0]; 


% non symmetric proportional constraint model

C{end+1} = [1,2,2,3,4,-1]; d = [d;0];
C{end+1} = [2,1,2,4,3,-1]; d = [d;0];
C{end+1} = [1,2,3,5,6,-1]; d = [d;0];
C{end+1} = [2,1,3,6,5,-1]; d = [d;0];

%{
% symmetry constraints for 20 modes
C{end+1} = [1,2,1,2,1,1]; d = 0;
C{end+1} = [3,4,1,4,3,1]; d = [d;0]; 
C{end+1} = [5,6,1,6,5,1]; d = [d;0];

C{end+1} = [1,2,2,3,4,-1]; d = [d;0];
C{end+1} = [1,2,3,5,6,-1]; d = [d;0];
%}
%{
C{end+1} = [7,8,1,8,7,1]; d = [d;0];
C{end+1} = [9,10,1,10,9,1]; d = [d;0];
C{end+1} = [11,12,1,12,11,1]; d = [d;0];
C{end+1} = [13,14,1,14,13,1]; d = [d;0];
C{end+1} = [15,16,1,16,15,1]; d = [d;0];
C{end+1} = [17,18,1,18,17,1]; d = [d;0];
C{end+1} = [19,20,1,20,19,1]; d = [d;0];
%}
%{

C{end+1} = [1,2,4,7,8,-1]; d = [d;0];
C{end+1} = [1,2,5,9,10,-1]; d = [d;0];
C{end+1} = [1,2,6,12,11,-1]; d = [d;0];
C{end+1} = [1,2,7,13,14,-1]; d = [d;0];
C{end+1} = [1,2,8,15,16,-1]; d = [d;0];
C{end+1} = [1,2,9,17,18,-1]; d = [d;0];
C{end+1} = [1,2,10,20,19,-1]; d = [d;0];
%}


%find best fit coefficients
if symmetry == 1
    %definde used values corrosponding to 
    range = [2:(size(X,2)/2-1) , (size(X,2)/2 +2) : (size(X,2)-1)];
    
    %use constrained sparsifying function

    %Xi = sparsifyDynamics(Theta(range,:),da,lambda,nVars);
    
    Xi = sparsifyDynamics_con(Theta(range,:),da,lambda2,nVars,C,d);

    %Xi = sparsifyDynamics_con_KKT(Theta(range,:),da,lambda2,nVars,C,d);
   
    %Xi = sparsifyDynamics_con_zero(Theta(range,:),da,lambda2,nVars,C);

    %Xi = sparsifyDynamics_con_single(Theta(range,:),da, lambda2, nVars,C,d);

    %Xi = sparsifyDynamics_con_lasso(Theta(range,:),da,lambda,nVars,C,d,alpha);

    %Xi = sparsifyDynamics_con_mix(Theta(range,:),da,lambda,lambda2,nVars,C,d,alpha);

    %Xi = sparsifyDynamics_con_mix2(Theta(range,:),da,lambda,lambda2,nVars,C);

else
    Xi = sparsifyDynamics(Theta(2:end-1,:),da,lambda,nVars);

    %Xi = sparsifyDynamics_con(Theta(2:end-1,:),da,lambda2,nVars,C,d);
    
    %Xi = sparsifyDynamics_con_mix(Theta(2:end-1,:),da,lambda,lambda2,nVars,C,d,alpha);

end


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
figure('position',[400,100,800,500])
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
    ylabel(sprintf('$a_{%i}$',i));
    
    if i == 1; legend('Reproduced','Original');
    elseif i == m; xlabel('Time [s]');
    end
    set(gca,'fontsize',9)
end

%plot phase space with discovered dynamics
pha_fig = figure('position',[400,100,800,500]); 
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
div_fig_name = "Circle Phase portrait.pdf"; 
exportgraphics(pha_fig, fullfile(savepath,div_fig_name), "ContentType","vector");

%% Recreate flow from identified system and POD modes 
savepath = "C:\Users\selia\Desktop\Bachelor Project images\Contraint analysis\Brunton data";
figtitle = "Original"; 
savename = append(figtitle," Snapshot.pdf"); 


%Sum all modes with discovered amplitudes. 
Recreate(:,:) = Xavg.*ones(size(X,1),length(tspan)); 
for j = 1:size(Recreate,2)
    for i = 1:m 
    Recreate(:,j) = Recreate(:,j) + U(:,i)*ai(j,i);
    end
end

timestep = 151; %chosen timestep for plotting

fig_recreated = figure('Position',[400,10,600,200]);
    %subplot(1,2,1);
    plotCylinder(reshape(X(:,timestep),nx,ny)); % plot of wake.
    set(gca,'fontsize',14); %title('Original'); 
    
    %subplot(1,2,2); 
    %plotCylinder(reshape(Recreate(:,timestep),nx,ny)); % plot of transformed wake.
    %set(gca,'fontsize',14); %title(figtitle);

exportgraphics(fig_recreated,fullfile(savepath,savename),"ContentType","vector")
%AnimateSquare(VORT_stack,'Square_flow_recreated_anim.gif',0.05,1);
%title('Recreated vorticity Snapshot')

%% plot enstropy with time

enstrophy_org = sum(X.^2,1); 
enstrophy_rec = sum(Recreate.^2,1); 

figure(); hold on; grid on; 
plot(t,enstrophy_rec); 
plot(t(Plt_inx),enstrophy_org(1:end/2));

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
