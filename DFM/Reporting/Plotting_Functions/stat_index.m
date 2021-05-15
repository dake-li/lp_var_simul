function index = stat_index(stat_name, settings)
% Function for giving the index of a particular summary statistic in all
% result matrices

% stat_name: 'mean', 'std', 'winsorized_mean', 'winsorized_std', 'quant_0.1', etc.
%         or simply input numbers like 0.1 for quantiles

% unpack
summ_stat_name = settings.simul.summ_stat_name;

% check if the stat_name is a number
if isnumeric(stat_name)
    stat_name = ['quant_',num2str(stat_name)]; % convert to a string
end

% find the index
index = find(strcmp(summ_stat_name, stat_name));

end

