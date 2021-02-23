function plot_choice(choice, cmap, horzs, weight_grid, methods_select, plot_name, plot_legend, legend_type, font_size)
% Function for plotting the optimal choice of method

figure
imagesc(horzs,weight_grid,choice)
colormap(cmap)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',7/8*font_size)
xlim([min(horzs) max(horzs)])
ylim([0 1])
set(gca,'XTick',horzs(mod(horzs,2) == 0));
set(gca,'YTick',[0:0.2:1]);
set(gca,'ydir','normal')
set(gca,'TickLength',[0 0])
xlabel('Horizon','interpreter','latex','FontSize',font_size);
ylabel('Bias Weight','interpreter','latex','FontSize',font_size);

if legend_type==0
    
    cbh = colorbar();
    caxis([1,max(methods_select)])
    set(cbh, 'YTick', 1:max(methods_select), ...
        'YTickLabel', plot_legend,'TickLabelInterpreter','latex')
    
else

    % Legend: trick from
    % https://www.mathworks.com/matlabcentral/answers/326478-add-legend-to-imagesc

    hidden_h = [];
    hold on;
    for K = 1:length(plot_legend)
        hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor', 'none');
    end
    hold off;
    colormap(cmap(1:length(plot_legend),:));
    uistack(hidden_h, 'bottom');
    legend(hidden_h, plot_legend, 'Location', 'eastoutside', 'NumColumns', 1, 'interpreter', 'latex','FontSize',font_size);

end

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) 1.4*pos(3) 1*pos(4)]);
set(gcf, 'PaperPositionMode', 'auto');

end