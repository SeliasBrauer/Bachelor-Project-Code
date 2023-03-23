%% Plot mode energies for multiple simulations. 
clc; clear variables; close all; 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% data from benajamin: 
ME_CON_LY4 = [23.6747322196424, 44.7409028101610, 53.2151766374616, 61.2099102127126 ...
    , 69.0434446518334, 76.8297184476719, 80.3988365283100, 83.9311659921872 ...
    , 87.1217782152595, 90.2618891044756];


ME_CON_LY6 = [24.7928240442085, 48.6605433933086, 58.1284990055368, ...
    67.4331862871242, 75.2428363802387, 82.8793891023216, 85.8854741418080 ...
    , 88.8840910008077, 91.2749270815915, 93.6338648220949];

%load ovn data: 
load Mode_energy.mat 

figure(); grid on; hold on; 

plot(ME_CYL,'*'); 
plot(ME_circ_vel,'*');
plot(ME_SQ_UNCON,'*');
plot(ME_CON_LY6,'*');  
plot(ME_CON_LY5,'*'); 
plot(ME_CON_LY4,'*'); 



 

legend('Cylinder vorticity','Cylinder velocity', 'Square unconfined','Square confined $L_y = 6$'... 
    ,'Square confined $L_y = 5$','Square confined $L_y = 4$',Location='southeast')

title('Energy caputed by the first 10 modes')
xlabel('Number of modes included')
ylabel('Energy captured by the model (\%)')
xlim([0 10]); ylim([0 100]); 