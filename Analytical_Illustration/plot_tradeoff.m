%% ANALYTICAL ILLUSTRATION: BIAS-VARIANCE TRADE-OFF
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clc
clear all
close all

addpath('./Subroutines');

%% SETTINGS

rhos = [0.6 0.9];   % Values of rho to consider
alphas = [1 5];     % Values of alpha to consider
sigma_2 = 1;        % Single fixed value for sigma_2
hs = 0:19;          % Horizons to plot

%% COMPUTE BIAS AND VARIANCE FOR LP AND VAR

for i_rho = 1:length(rhos)
    rho = rhos(i_rho);
    alpha = alphas(i_rho);
    [bias_var(i_rho,:), var_var(i_rho,:), var_lp(i_rho,:)] = asy_bias_var(rho',sigma_2,alpha,hs);
end

bias_lp = 0 * bias_var;

%% PLOT RESULTS

% settings

line_colors = lines(7);
lp_indx  = 2; % LP color
var_indx = 4; % VAR color
linespecs = {'-', '-o'};
line_width = 3;
marker_size = 3.5;

plotwidth = 0.4;
gapsize = 0.05;
gapsize_edges = (1-2*plotwidth-1*gapsize)/2;
left_pos = [gapsize_edges, gapsize_edges + gapsize + plotwidth];

% bias/variance figure

figure(1);

subplot(1,2,1)
pos = get(gca, 'Position');
pos(1) = left_pos(1);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'TickLabelInterpreter','latex')
for i_rho = 1:length(rhos)
    hold on
    plot(hs,bias_lp(i_rho,:), linespecs{i_rho}, 'Color', line_colors(lp_indx,:), 'LineWidth', line_width, 'MarkerSize', marker_size);
end
for i_rho = 1:length(rhos)
    hold on
    plot(hs,bias_var(i_rho,:), linespecs{i_rho}, 'Color', line_colors(var_indx,:), 'LineWidth', line_width, 'MarkerSize', marker_size);
end
set(gcf,'color','w')
xlim([min(hs) max(hs)]);
xticks(min(hs):4:max(hs))
ylim([0 3.5])
title('Bias','interpreter','latex')
xlabel('Horizon','interpreter','latex')
legend({'LP','LP','VAR $\,\quad \{ \rho = 0.6,\, \alpha = 1\}$', ...
    'VAR $\,\quad \{ \rho = 0.9,\, \alpha = 5\}$'},'Location','North','NumColumns', 2, 'interpreter','latex')
grid on
hold off
set(gca, 'FontSize', 13);
set(gca, 'TitleFontSizeMultiplier', 1.2);

subplot(1,2,2)
pos = get(gca, 'Position');
pos(1) = left_pos(2);
pos(3) = plotwidth;
set(gca,'Position', pos)
set(gca,'TickLabelInterpreter','latex')
for i_rho = 1:length(rhos)
    hold on
    plot(hs,sqrt(var_lp(i_rho,:)), linespecs{i_rho}, 'Color', line_colors(lp_indx,:), 'LineWidth', line_width, 'MarkerSize', marker_size);
end
for i_rho = 1:length(rhos)
    hold on
    plot(hs,sqrt(var_var(i_rho,:)), linespecs{i_rho}, 'Color', line_colors(var_indx,:), 'LineWidth', line_width, 'MarkerSize', marker_size);
end
set(gcf,'color','w')
xlim([min(hs) max(hs)]);
xticks(min(hs):4:max(hs))
ylim([0 3.5])
title('Standard Deviation','interpreter','latex')
xlabel('Horizon','interpreter','latex')
grid on
hold off
set(gca, 'FontSize', 13);
set(gca, 'TitleFontSizeMultiplier', 1.2);

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.8*0.85*pos(3) 1.8*0.425*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');
print('biasvar_simple','-depsc');