%----------------------------------------------------------------
% File/Folder Names
%----------------------------------------------------------------

exper_filename = exper_files{ne}; % Name of current experiment
exper_plotname = exper_names{ne};
file_name = fullfile(lags_folders{nf}, exper_filename); % Name of .mat results file
output_folder = fullfile(output_dir, file_name); % Name of output folder
mkdir(output_folder); % Create output folder        

%----------------------------------------------------------------
% Load Simulation Results
%----------------------------------------------------------------

% see if this exper is start of group
if (ne == 1) || (exper_group_end(ne-1) == 1)
    res = load(fullfile(rootfolder, strcat(file_name, '.mat'))); % Directly load
else
    res_part = load(fullfile(rootfolder, strcat(file_name, '.mat'))); % Load
    res = combine_exper(res, res_part); % Merge
end

horzs = res.settings.est.IRF_select; % Impulse response horizons
methods_names_plot = methods_names{ne};