function plot_tradeoff(pref_base, cmap, horzs, weight_grid, plot_name, font_size)
% Function for the trade-off between two methods across DGP

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
xlabel('Horizon','interpreter','latex','FontSize',font_size);
ylabel('Bias Weight','interpreter','latex','FontSize',font_size);

cbh = colorbar();
caxis([0,1])
set(cbh, 'YTick', [0:0.2:1], ...
    'YTickLabel', [0:0.2:1],'TickLabelInterpreter','latex')

end