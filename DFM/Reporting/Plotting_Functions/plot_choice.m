function plot_choice(choice, cmap, horzs, weight_grid, methods_select, plot_name, plot_legend, legend_type)

figure
imagesc(horzs,weight_grid,choice)
colormap(cmap)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
xlim([min(horzs) max(horzs)])
ylim([0 1])
set(gca,'XTick',horzs(mod(horzs,2) == 0));
set(gca,'YTick',[0:0.2:1]);
set(gca,'ydir','normal')
set(gca,'TickLength',[0 0])
title(plot_name, 'interpreter', 'latex','FontSize',20);

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
    legend(hidden_h, plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'interpreter', 'latex','FontSize',16);

end

end