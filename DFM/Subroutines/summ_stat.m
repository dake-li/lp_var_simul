function stat_vector = summ_stat(x_vector, winsor_percent, quantiles)
% summarize the dim of MC into a few stats: mean/std, winsorized mean/std,
% and quantiles. Return the stat_vector in the similar shape of x_vector

size_vector = size(x_vector); % record the dim of MC
n_MC = length(x_vector);
x_vector = squeeze(x_vector); % n_MC by 1 vector

this_mean = mean(x_vector);
this_std = std(x_vector);
this_winsorized_x_vector = rmoutliers(x_vector, 'percent', [100*winsor_percent, 100*(1-winsor_percent)]);
this_winsorized_mean = mean(this_winsorized_x_vector);
this_winsorized_std = std(this_winsorized_x_vector);
this_quantiles = quantile(x_vector, quantiles);

stat_vector = [this_mean, this_std, this_winsorized_mean, this_winsorized_std, this_quantiles];
size_vector(size_vector==n_MC) = length(stat_vector);
stat_vector = reshape(stat_vector, size_vector);

end

