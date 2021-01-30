function add_legend(ax_vector, plot_legend, font_size)

    % Function for adding legend
    
    for i = 1:length(ax_vector)
        ax_vector(i).Position(2) = 0.25;
        ax_vector(i).Position(4) = 0.6;
    end
    legend(plot_legend, 'Position', [0.45,0.05,0.1,0.05], 'NumColumns', 3, 'interpreter', 'latex','FontSize',font_size);

end
