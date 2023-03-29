clc; clear variables; close all; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%Import data
vorticity = readmatrix('C:\Users\selia\Downloads\vorticity-1D11D4D.csv'); 
Lift_force = readmatrix('C:\Users\selia\Downloads\F_y-1D11D4D.csv');


% cutoff data before settling
vorticity(1:300,:) = []; 
Lift_force(1:300,:) = []; 

%plot data:
figure()
plot(vorticity(:,1),vorticity(:,2))

figure()
plot(Lift_force(:,1),Lift_force(:,2))

%define timestep; 
dt = Lift_force(2,1)-Lift_force(1,1);

%strings for titles (Optional) 
string1 = 'Lift force frequency sprectrum';
string2 = 'Vorticity frequency sprectrum';

%return largest frequency responce
st1 = Strouhal_number_f(Lift_force(:,2),dt)
st2 = Strouhal_number_f(vorticity(:,2),dt);


function st = Strouhal_number_f(Data, dt)
% domain properties; 
D = 1; U = 1; 
[pks,locs] = findpeaks(Data); %find location of peaks

%find # timesteps between peaks
timestep_period = mean(locs(2:end) - locs(1:end-1));

st = 1/(timestep_period*dt) * D/U;
end





%{
function st = Strouhal_number_f(Data,dt,plottitle)
%calculates the strouhal number from collected data by finding the largest 
%frequency responce by an FFT. 
%Takes data vector, timestep and title string for plotting. 

D = 1; U = 1; %Flow domain properties

figure()
T = length(Data)*dt;  % length of signal
fs = 1/dt; %sample frequency

Fmax = 0.5*fs; %max dectecable frequency
Fmin = 1/T; %min detectable frequency

Y = fft(Data) / length(Data); %the actual fft

Y(1) = []; % for plotting only
Yplot = Y(1: floor(length(Data)/2)); 
Power = abs(Yplot).^2; 
fspan = linspace(Fmin,Fmax, floor(length(Data)/2)); % frequency span

bar(fspan,Power);
xlabel('Freqency (Hz)'); ylabel('Power')
title(plottitle); %makes title of plot. 

[~,I_max] = max(Power); % find index of maximum responce
dom_freq = fspan(I_max); % find dominant frequency 
st = dom_freq*D / U; 

end
%}
