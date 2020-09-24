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

res = load(fullfile(rootfolder, strcat(file_name, '.mat'))); % Load
horzs = res.settings.est.IRF_select; % Impulse response horizons
methods_names_plot = methods_names{ne};