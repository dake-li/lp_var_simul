function plot_tradeoff(pref_base, cmap, horzs, weight_grid, plot_name, font_size)

figure
imagesc(horzs,weight_grid,pref_base)
colormap(cmap)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',7/8*font_size)
xlim([min(horzs) max(horzs)])
ylim([0 1])
set(gca,'XTick',horzs(mod(horzs,2) == 0));
set(gca,'YTick',[0:0.2:1]);
set(gca,'ydir','normal')
set(gca,'TickLength',[0 0])
title(plot_name, 'interpreter', 'latex','FontSize',9/8*font_size);
xlabel('Horizon','interpreter','latex','FontSize',font_size);
ylabel('Bias Weight','interpreter','latex','FontSize',font_size);

cbh = colorbar();
caxis([0,1])
set(cbh, 'YTick', [0:0.2:1], ...
    'YTickLabel', [0:0.2:1],'TickLabelInterpreter','latex')

% figure
% % h = heatmap(pref_base,'ColorbarVisible','off','GridVisible','off');
% h = heatmap(pref_base,'GridVisible','off');
% h.Colormap = cmap;
% h.ColorLimits = [0 1];
% h.CellLabelColor = 'none';
% CustomXLabels = string(horzs);
% CustomXLabels(:) = " ";
% h.XDisplayLabels = CustomXLabels';
% CustomYLabels = weight_grid;
% CustomYLabels(:) = " ";
% h.YDisplayLabels = CustomYLabels';
% set(struct(h).NodeChildren(3), 'XTickLabelRotation', 0);
% 
% a2 = axes('Position', h.Position);
% set(gca,'TickLabelInterpreter','latex')
% a2.FontSize = 14;
% xlim([min(horzs) max(horzs)])
% ylim([0 1])
% a2.Color = 'none';  
% a2.Title = title(plot_name, 'interpreter', 'latex','FontSize',20);
% xlabel('Horizon','interpreter','latex','FontSize',16);
% ylabel('Bias Weight','interpreter','latex','FontSize',16);
% a2.XTick = horzs(mod(horzs,2) == 0);
% a2.YTick = [0:0.2:1];

end