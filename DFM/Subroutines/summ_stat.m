function stat_vector = summ_stat(x_vector, winsor_percent, quantiles)
% Function for summarizing the result matrices
    % summarize the dim of MC into a few stats:
        % mean/std
        % winsorized mean/std
        % quantiles
    
    % return the stat_vector in the similar shape of x_vector

size_vector = size(x_vector); % record the dim of MC
n_MC = length(x_vector);
x_vector = squeeze(x_vector); % n_MC by 1 vector

this_mean = mean(x_vector);
this_std = std(x_vector);

this_quantiles = quantile(x_vector, quantiles);

lower_bound = quantile(x_vector, winsor_percent);
upper_bound = quantile(x_vector, 1 - winsor_percent);
this_winsorized_x_vector = max(min(x_vector, upper_bound), lower_bound);
this_winsorized_mean = mean(this_winsorized_x_vector);
this_winsorized_std = std(this_winsorized_x_vector);


stat_vector = [this_mean, this_std, this_winsorized_mean, this_winsorized_std, this_quantiles];
size_vector(size_vector==n_MC) = length(stat_vector);
stat_vector = reshape(stat_vector, size_vector);

end

