function plot_frontier(bias, std, bias_frontier_grid, max_std, method_labels, font_size, plot_title)
% Function for plotting the bias/variance frontier

    % settings

    line_colors = repmat(lines(7),2,1);
    colors_indx = [4 3 1 2 5 6 7];
    line_colors = line_colors(colors_indx,:);
    
    % compute lower envelope for visual aid
    
    std_frontier_grid = std_frontier_fn(bias,std,bias_frontier_grid);
    
    % figures
    
    set(gca,'FontSize',7/8*font_size)
    set(gca,'TickLabelInterpreter','latex')
    
    grid on;
    
    % Points
    for i = 1:length(bias)
        hold on;
        plot(bias(i), std(i), '.', 'Color', line_colors(i,:), 'MarkerSize', 30);
    end
    
    % Frontier curve for visual aid
    hold on;
    plot(bias_frontier_grid,std_frontier_grid,'Color',[0.4 0.4 0.4],'LineWidth',2,'LineStyle','--')
    hold off
    
    xlabel('Abs. Bias (Median Across DGPs)','interpreter','latex','FontSize',font_size);
    ylabel('Std. Dev. (Median Across DGPs)','interpreter','latex','FontSize',font_size);
    title(plot_title,'interpreter','latex','FontSize',1.1*font_size);
    xlim([0 max(bias_frontier_grid)]);
    ylim([0 max_std]);
    
    % Text annotations
    [xnorm, ynorm] = coord2norm(gca, bias, std);
    for i = 1:length(bias)
        annotation('textarrow', xnorm(i)+[method_labels{i,2}(1) 0], ynorm(i)+[method_labels{i,2}(2) 0], ...
                   'String', method_labels{i,1}, 'Interpreter', 'latex', 'FontSize', 0.9*font_size, 'TextMargin', 1, ...
                   'TextBackgroundColor', 'w', 'TextEdgeColor', [0.6 0.6 0.6], 'HeadLength', 5, 'HeadWidth', 5);
    end

end