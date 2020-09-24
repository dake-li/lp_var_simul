function plot_choice(choice, horzs, weight_grid, plot_name, plot_legend)

figure
imagesc(horzs,weight_grid,choice)
colormap(lines)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',14)
xlim([min(horzs) max(horzs)])
ylim([0 1])
set(gca,'XTick',horzs(mod(horzs,2) == 0));
set(gca,'YTick',[0:0.2:1]);
set(gca,'ydir','normal')
set(gca,'TickLength',[0 0])
title(plot_name, 'interpreter', 'latex','FontSize',20);


% Legend: trick from
% https://www.mathworks.com/matlabcentral/answers/326478-add-legend-to-imagesc

hidden_h = [];
hold on;
for K = 1:length(plot_legend)
    hidden_h(K) = surf(uint8([K K;K K]), 'edgecolor', 'none');
end
hold off;
colormap(lines(length(plot_legend)));
uistack(hidden_h, 'bottom');
legend(hidden_h, plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'interpreter', 'latex','FontSize',16);

end