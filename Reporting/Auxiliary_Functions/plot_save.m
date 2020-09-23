function plot_save(filename, suffix)

    % Function for saving figures
    
    saveas(gcf, strcat(filename, '.', suffix));
    close(gcf);

end