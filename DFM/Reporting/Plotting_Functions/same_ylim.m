function same_ylim(ax_vector)
    % Enforce same ylim
    if length(ax_vector) == 1
        return;
    end
    
    ylim = NaN(length(ax_vector),2);
    for i = 1:length(ax_vector)
        ylim(i,:) = get(ax_vector(i), 'YLim');
    end
    ylim_new = [min(ylim(:,1)) max(ylim(:,2))];
    for i = 1:length(ax_vector)
        set(ax_vector(i), 'YLim', ylim_new);
    end

end
