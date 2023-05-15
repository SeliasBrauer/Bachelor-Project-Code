clc; clear variables; close all; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

load Derivative_error.mat
dt = [0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2]; 

dt = 1./dt;

error(1,:) = div_mse_001;
error(2,:) = div_mse_002;
error(3,:) = div_mse_003;
error(4,:) = div_mse_005;
error(5,:) = div_mse_008;
error(6,:) = div_mse_01;
error(7,:) = div_mse_015;
error(8,:) = div_mse_02;

t = ones(1,8); 
figure; hold on; grid on; 

for i = 1:6
 scatter(dt, error(:, i),'filled',DisplayName=sprintf('Mode %i',i))
end 
legend('show')
xlabel('Sampling rate [Hz]'); 
ylabel('Mean squared error');
title('Derivative error (confined)')

%% convergance of Finite difference confined: 
clc; close all;  clear variables
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

dt = [0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.1, 0.15, 0.2]; 
dt = 1./dt;
load Derivatives_confined.mat % load data and add to large struct
data(1) = struct_001; 
data(2) = struct_002; 
data(3) = struct_003; 
data(4) = struct_004; 
data(5) = struct_005; 
data(6) = struct_008; 
data(7) = struct_010; 
data(8) = struct_015; 
data(9) = struct_020; 


% flip derivatives :
data(2).derivative(:,6) = -data(2).derivative(:,6) ;

data(3).derivative(:,4) = -data(3).derivative(:,4) ;
data(3).derivative(:,5) = -data(3).derivative(:,5) ;

data(4).derivative(:,6) = -data(4).derivative(:,6) ;

data(5).derivative(:,3) = -data(5).derivative(:,3) ;
data(5).derivative(:,4) = -data(5).derivative(:,4) ;
data(5).derivative(:,5) = -data(5).derivative(:,5) ;

data(6).derivative(:,3) = -data(6).derivative(:,3) ;
data(6).derivative(:,5) = -data(6).derivative(:,5) ;

data(7).derivative(:,3) = -data(7).derivative(:,3) ;
data(7).derivative(:,6) = -data(7).derivative(:,6) ;

data(8).derivative(:,3) = -data(8).derivative(:,3) ;
data(8).derivative(:,6) = -data(8).derivative(:,6) ;

data(9).derivative(:,3) = -data(9).derivative(:,3) ;
data(9).derivative(:,6) = -data(9).derivative(:,6) ;


% calculate error; 
figure(); hold on; 
for i = 1:length(data) 
%find intersecting timesteps
[C,ia,ib] = intersect(round(data(1).timeseries(2:end-1),2), round(data(i).timeseries(2:end-1),2));

% load derivatives: 
dt_true = data(1).derivative(1:end/2,:);
dt_esti= data(i).derivative(1:end/2,:);

%{
figure(); 
plot(dt_true(ia,6)); hold on 
plot(dt_esti(ib,6)); hold off
%}

mse(i,:) = 1/size(dt_true(ia,:),1) .* sum ( (dt_esti(ib,:) - dt_true(ia, :)).^2 ,1);
mse2(i,:) = rmse(dt_esti(ib,:),dt_true(ia,:));

difference = dt_esti(ib,:) - dt_true(ia, :); 

subplot(5,2,i)
plot(data(i).timeseries(2:end-1),difference)

if i == 1
    legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6')
end
xlabel('Time [s]'); ylabel('Difference')

p_div(i,:) =   abs(mean( (dt_esti(ib,:) - dt_true(ia, :)) ./ dt_true(ia,:) *100 ,1));

end


figure; hold on; grid on; 

for i = 1:6
 plot(dt, mse2(:, i),'--o',DisplayName=sprintf('Mode %i',i))
end 
set(gca,'yscale','log')
xlim([0,51])
legend('show')
xlabel('Sampling rate [Hz]'); 
ylabel('Root Mean Square Error');
title('Derivative error (confined)')


%% convergance of Finite difference unconfined: 
clc; close all;  clear variables
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

dt = 0.02 : 0.02 : 0.20;
dt = 1./dt;
load Derivatives_unconfined.mat % load data and add to large struct
data(1) = struct_002; 
data(2) = struct_004; 
data(3) = struct_006; 
data(4) = struct_008; 
data(5) = struct_010; 
data(6) = struct_012; 
data(7) = struct_014; 
data(8) = struct_016; 
data(9) = struct_018;
data(10) = struct_020;


% flip derivatives :
data(2).derivative(:,5) = -data(2).derivative(:,5) ;
data(2).derivative(:,6) = -data(2).derivative(:,6) ;

data(3).derivative(:,3) = -data(3).derivative(:,3) ;
data(3).derivative(:,5) = -data(3).derivative(:,5) ;

data(4).derivative(:,5) = -data(4).derivative(:,5) ;

data(5).derivative(:,4) = -data(5).derivative(:,4) ;
data(5).derivative(:,5) = -data(5).derivative(:,5) ;

data(6).derivative(:,3) = -data(6).derivative(:,3) ;
data(6).derivative(:,5) = -data(6).derivative(:,5) ;

data(9).derivative(:,5) = -data(9).derivative(:,5) ;

data(10).derivative(:,5) = -data(10).derivative(:,5) ;


% calculate error; 
figure(); hold on; 
for i = 1:length(data) 
%find intersecting timesteps
[C,ia,ib] = intersect(round(data(1).timeseries(2:end-1),2), round(data(i).timeseries(2:end-1),2));

% load derivatives: 
dt_true = data(1).derivative(1:end/2,:);
dt_esti= data(i).derivative(1:end/2,:);

%{
figure(); 
plot(dt_true(ia,6)); hold on 
plot(dt_esti(ib,6)); hold off
%}

mse(i,:) = 1/size(dt_true(ia,:),1) .* sum ( (dt_esti(ib,:) - dt_true(ia, :)).^2 ,1);
mse2(i,:) = rmse(dt_esti(ib,:),dt_true(ia,:));

difference = dt_esti(ib,:) - dt_true(ia, :); 

subplot(5,2,i)
plot(data(i).timeseries(2:end-1),difference)

if i == 1
    legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6')
end
xlabel('Time [s]'); ylabel('Difference')

p_div(i,:) =   abs(mean( (dt_esti(ib,:) - dt_true(ia, :)) ./ dt_true(ia,:) *100 ,1));

end


figure; hold on; grid on; 

for i = 1:6
 plot(dt, p_div(:, i),'--o',DisplayName=sprintf('Mode %i',i))
end 
set(gca,'yscale','log')
xlim([0,51])
legend('show')
xlabel('Sampling rate [Hz]'); 
ylabel('Percent deviation (\%)');
title('Derivative error (unconfined)')





