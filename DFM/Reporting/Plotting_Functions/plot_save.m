function plot_save(filename, suffix)

    % Function for saving figures

    if ~iscell(suffix)
        suffix = {suffix};
    end

    for i = 1:length(suffix) % go across multiple suffixes

        the_suffix = suffix{i};

        if the_suffix == 'eps'
            saveas(gcf, strcat(filename, '.', the_suffix), 'epsc');
        else
            saveas(gcf, strcat(filename, '.', the_suffix));
        end

    end
    
    close(gcf);

end