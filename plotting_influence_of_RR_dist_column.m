% plotting_influence_of_RR_dist_column
clc; clear all; close all; 
[~, ~, raw] = xlsread('/Users/Clemens/Documents/Uni/6. Semester/Case Study II/Influence of changing RR_multiplier.xlsx','Tabelle1');
raw = raw(3:12,:);
data = reshape([raw{:}],size(raw));
RR = data(:,1);
height = data(:,2);
CAPEX = data(:,3);
OPEX = data(:,4);
radius = data(:,5);
clearvars data raw;
figure;
subplot(2, 2, 1);  
fontsize1 = 14; 
plot(RR, height, 'Color', 'Blue', 'LineWidth', 1.5);
title('Column height vs. Reflux ratio', 'Fontsize', fontsize1);
xlabel('$RR factor$', 'Fontsize', fontsize1);
ylabel('Column height $H$ [m]', 'Fontsize', fontsize1);

subplot(2,2, 2);       
plot(RR, 2*radius, 'Color', 'Blue', 'LineWidth', 1.5);     
title('Diameter vs. Reflux ratio', 'Fontsize', fontsize1);
xlabel('$RR factor$', 'Fontsize', fontsize1);
ylabel('Column diameter $d$ [m]', 'Fontsize', fontsize1);

subplot(2,2, 3);       
plot(RR, CAPEX*1e-6, 'Color', 'Blue', 'LineWidth', 1.5);      
title('CAPEX for Column vs. Reflux ratio', 'Fontsize', fontsize1);
xlabel('$RR factor$', 'Fontsize', fontsize1);
ylabel('CAPEX [M\$]', 'Fontsize', fontsize1);

subplot(2,2, 4);       
plot(RR, OPEX*1e-6, 'Color', 'Blue', 'LineWidth', 1.5);      
title('OPEX for Column vs. Reflux ratio', 'Fontsize', fontsize1);
xlabel('$RR factor$', 'Fontsize', fontsize1);
ylabel('OPEX [M\$/a]', 'Fontsize', fontsize1);