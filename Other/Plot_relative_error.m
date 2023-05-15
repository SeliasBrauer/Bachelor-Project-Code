clc; clear variables; close all; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

load Error_struct.mat

data(1) = error_struct_hierarchical; 
data(2) = error_struct_hierarchical2; 
data(3) = error_struct_pair;
data(4) = error_struct_pair2; 
data(5) = error_struct_symmetric; 
data(6) = error_struct_symmetric2; 


figure(); hold on; grid on; 
for i = 1:length(data)
    string = append(data(i).string, sprintf('~~dt = %.2f s',data(i).dt)); 
    plot(data(i).re,'--o','DisplayName',string); 
end 
legend show; 
legend Location northwest; 
xlabel('Mode number'); ylabel('Amplitude relative RMSE'); 

%{
axes('position',[0.16 .65 .3 .3]); hold on; grid on; 
box on % put box around new pair of axes
for i = 1:length(data)
    indexOfInterest = 1:6; % range of t near perturbation
    plot(1:6,data(i).re(indexOfInterest),'--o') % plot on new axes
    axis tight

end 
%}