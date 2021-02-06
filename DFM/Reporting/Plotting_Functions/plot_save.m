function plot_save(filename, suffix)

    % Function for saving figures
    
    if suffix == 'eps'
        saveas(gcf, strcat(filename, '.', suffix), 'epsc');
    else
        saveas(gcf, strcat(filename, '.', suffix));
    end
    close(gcf);

end