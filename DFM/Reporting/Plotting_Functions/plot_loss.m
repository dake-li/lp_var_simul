function plot_loss(horzs, results, add_line, plot_name, plot_legend, font_size, varargin)
% Function for plotting bias/variance loss results across estimation methods

    % settings

    line_colors = repmat(lines(7),2,1);
    colors_indx = [4 3 1 2 5 6 7];
    line_colors = line_colors(colors_indx,:);
    line_styles = {'-', '--', '-x', '-', '-.', '-o', ':'};
    line_width  = repmat(2,14,1);
    line_width(1) = 5;
    line_width(4) = 5;
    
    % figures
    
    if isempty(varargin)
        figure;
    end
    set(gca,'FontSize',7/8*font_size)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    
    hold on;
    for i=1:size(results,2)
        plot(horzs, results(:,i)', line_styles{i}, 'Color', line_colors(i,:), 'LineWidth', line_width(i));
        xlim([min(horzs) max(horzs)])
    end
    hold off;
    
    if ~isempty(add_line)
        hold on;
        the_xlim = xlim;
        line(the_xlim, add_line*ones(1,2), 'Color', 'k', 'LineStyle', ':');
        hold off;
        xlim(the_xlim);
    end
    xlabel('Horizon','interpreter','latex','FontSize',font_size);
    set(gca, 'XTick', [min(horzs) 2:2:max(horzs)]);
    if isempty(varargin)
        legend(plot_legend, 'Location', 'eastoutside', 'NumColumns', 1, 'interpreter', 'latex','FontSize',font_size);
    end
    
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');

end