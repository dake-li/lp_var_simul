%% ANALYTICAL ILLUSTRATION: INDIFFERENCE WEIGHTS
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clc
clear all
close all

addpath('./Subroutines');

%% SETTINGS

rhos = [0.2 0.6 0.9];   % Values of rho to consider
alphas = [1 5 10];      % Values of alpha to consider
sigma_2 = 1;            % Single fixed value for sigma_2
hs = 2:19;              % Horizons to plot

linespecs = {'-', '-x', '-*'}; % Linespecs for different alpha values

%% INDIFFERENCE WEIGHT PLOT

% settings

num_horz = length(hs);
num_rho = length(rhos);
num_alpha = length(alphas);
the_legend = cellfun(@(x) sprintf('%s%2d%s', '$\alpha=', x, '$'), num2cell(alphas), 'UniformOutput', false); % Figure legend

line_colors = repmat([0 0.3 0.6]',1,3);

% computations and figure constructions

figure('Units', 'inches', 'Position', [0 0 10 3]);
for ir=1:num_rho
    
    subplot(1,num_rho,ir);
    set(gca,'TickLabelInterpreter','latex')
    
    % Compute asy. bias/var for different alphas and horizons
    [bias_var, var_var, var_lp] = asy_bias_var(rhos(ir),sigma_2,alphas',hs);
    
    % Loss function is w*bias^2 + (1-w)*var
    % Compute smallest weight w such that LP is preferred over VAR
    w_indiff = (var_lp-var_var)./(bias_var.^2+var_lp-var_var);
    
    % Plot indifference weights across horizons, separately for each alpha
    hold on;
    for ia=1:num_alpha
        plot(hs,w_indiff(ia,:),linespecs{ia},'LineWidth',2,'Color',line_colors(ia,:));
    end
    hold off;
    title(sprintf('%s%3.1f%s', '$\rho=', rhos(ir), '$'), 'Interpreter', 'latex');
    grid on;
    xlim([min(hs) max(hs)]);
    xlabel('Horizon','Interpreter','latex');
    xticks([min(hs) 4:4:max(hs)])
    ylim([0 1]);
    set(gca, 'FontSize', 12);
    set(gca, 'TitleFontSizeMultiplier', 1.2);
    
    % Add legend
    if ir==1
        legend(the_legend, 'Interpreter', 'latex', 'Location', 'SouthEast', 'FontSize', 12);
    end
    
end
print('weight_indiff','-depsc');