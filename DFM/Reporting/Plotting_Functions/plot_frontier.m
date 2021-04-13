function plot_frontier(loss, frontier, methods, font_size, plot_title, pos, dist)
% Function for plotting the bias/variance frontier

    % settings

    line_colors = repmat(lines(7),2,1);
    colors_indx = [4 3 1 2 5 6 7];
    line_colors = line_colors(colors_indx,:);
    
    % extract inputs
    
    bias_methods  = loss(:,1);
    std_methods   = loss(:,2);
    bias_frontier = frontier(:,1);
    std_frontier  = frontier(:,2);
    
    % figures
    
    figure
    set(gca,'FontSize',7/8*font_size)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    hold on
    for i = 1:size(bias_methods,1)
%         text(bias_methods(i),std_methods(i),methods{i},'Color', line_colors(i,:),'FontSize',20,'interpreter','latex','FontWeight','bold');
        plot(bias_methods(i),std_methods(i),'x','Color', line_colors(i,:),'MarkerSize',20,'LineWidth',10)
        if pos(i) == 1
        labelpoints(bias_methods(i),std_methods(i),methods{i},'NE',dist(i),'interpreter','latex','FontSize',font_size)
        elseif pos(i) == 2
        labelpoints(bias_methods(i),std_methods(i),methods{i},'SE',dist(i),'interpreter','latex','FontSize',font_size)
        elseif pos(i) == 3
        labelpoints(bias_methods(i),std_methods(i),methods{i},'SW',dist(i),'interpreter','latex','FontSize',font_size)
        elseif pos(i) == 4
        labelpoints(bias_methods(i),std_methods(i),methods{i},'NW',dist(i),'interpreter','latex','FontSize',font_size)
        end
        hold on
    end
    hold on
    plot(bias_frontier,std_frontier,'Color',[0.4 0.4 0.4],'LineWidth',4,'LineStyle','--')
    xlabel('Bias','interpreter','latex','FontSize',font_size);
    ylabel('St. Dev.','interpreter','latex','FontSize',font_size);
    title(plot_title,'interpreter','latex','FontSize',5/4*font_size);
%     legend(methods, 'Location', 'eastoutside', 'NumColumns', 1, 'interpreter', 'latex','FontSize',font_size);
    ylim([0 max(std_methods) * 1.1])
    hold off
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1.3*pos(4)]);
    set(gcf, 'PaperPositionMode', 'auto');

end