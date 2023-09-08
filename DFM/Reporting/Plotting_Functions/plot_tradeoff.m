function plot_tradeoff(pref_base, cmap, horzs, weight_grid, n_bin, font_size)
% Function for the trade-off between two methods across DGP

figure
imagesc(horzs,weight_grid,pref_base)
colormap(cmap)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',7/8*font_size)
xlim([min(horzs) max(horzs)])
ylim([0.5 1])
set(gca,'XTick',horzs(mod(horzs,2) == 0));
set(gca,'YTick',[0.5:0.1:1]);
set(gca,'ydir','normal')
set(gca,'TickLength',[0 0])
xtickangle(0);
xlabel('Horizon','interpreter','latex','FontSize',font_size);
ylabel('Bias Weight','interpreter','latex','FontSize',font_size);

% Legend
cbh = colorbar();
clim([0,1])
the_legend = cell(1,n_bin);
the_mult = 100/n_bin;
the_legend{1} = sprintf('$<%d%s$',the_mult,'%');
the_legend{end} = sprintf('$>%d%s$',the_mult*(n_bin-1),'%');
for j=2:n_bin-1
    the_legend{j} = sprintf('%d-%d%s',(j-1)*the_mult,j*the_mult,'%');
end
set(cbh, 'YTick', 0.5/n_bin:1/n_bin:1-0.5/n_bin, ...
    'YTickLabel', {'$<20$\%','20-40\%','40-60\%','60-80\%','$>80$\%'},'TickLabelInterpreter','latex','TickLength',0)

end