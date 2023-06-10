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
title('Unconfined case')

%{
axes('position',[0.16 .65 .3 .3]); hold on; grid on; 
box on % put box around new pair of axes
for i = 1:length(data)
    indexOfInterest = 1:6; % range of t near perturbation
    plot(1:6,data(i).re(indexOfInterest),'--o') % plot on new axes
    axis tight

end 
%}
%% Plot relative error obtained from data sensitivity analysis. 
load error_struct_data_analysis.mat
% figure save path
path = "C:\Users\selia\Desktop\Bachelor Project images\Data sensitivity analysis\Comparison plots";
f1_name = "Sampling rate dependence.pdf";
f2_name = "Total time dependence.pdf";

%plot singular value Relative RMSE 
f1 = figure(); hold on; grid on; 
f1.Units = 'centimeters'; 
f1.Position = [1,1,17,8];
hold on; grid on; 
order = [16,1,4,7,17,2,5,8,18,3,6,9];
symbol = {"--ko","--ksquare","--k^","-.ko","-.k square","-.k^",":ko",":ksquare",":k^",".-ko",".-ksquare",".-k^","--ko","--ksquare","--k^","-ko","-k square","-k^",":ko",":ksquare",":k^",".-ko",".-ksquare",".-k^"};

for i = order
    string = append(error_struct_data_analysis.string{i}, sprintf('~~dt = %.2f s',error_struct_data_analysis.dt(i))); 
     
    plot(error_struct_data_analysis.re{i},symbol{i},'DisplayName',string);
end
legend show; 
lgd = legend;
lgd.Location = 'northeastoutside'; 
xticks([1,2,3,4,5,6]);
set(gcf,"Color",'w')
xlabel('Mode number'); ylabel('Normalized RMSE');


exportgraphics(f1,fullfile(path, f1_name),'ContentType','vector')



%plot singular value Relative RMSE 
f2 = figure(); hold on; grid on; 
f2.Units = 'centimeters'; 
f2.Position = [1,1,17,8];
order = [13,10,16,14,11,17,15,12,18];
symbol = {"--ko","--ksquare","--k^","-ko","-k square","-k^",":ko",":ksquare",":k^"};

for i = order
    string = append(error_struct_data_analysis.string{i}, sprintf('~~T = %.2f s',error_struct_data_analysis.T(i))); 
    plot(error_struct_data_analysis.re{i},symbol{i-9},'DisplayName',string);
end
legend show; 
lgd = legend;
lgd.Location = 'northeastoutside'; 
xticks([1,2,3,4,5,6]);
set(gcf,"Color",'w')
xlabel('Mode number'); ylabel('Normalized RMSE');

exportgraphics(f2,fullfile(path, f2_name),'ContentType','vector')


