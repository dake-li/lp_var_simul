function plot_choice(choice, cmap, horzs, weight_grid, methods_select, plot_name)

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

cbh = colorbar();
caxis([1,max(methods_select)])
set(cbh, 'YTick', [1 2 3 4 5 6 7], ...
    'YTickLabel', {'LP', 'Pen. LP', 'BC-VAR','VAR','VAR-avg.','BVAR','SVAR-IV'},'TickLabelInterpreter','latex')

end