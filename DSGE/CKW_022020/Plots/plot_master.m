%% LP vs VAR: DFM SIMULATION STUDY
% Dake Li and Christian Wolf
% this version: 12/12/2019

%% HOUSEKEEPING

clc
clear all
close all

% path = 'D:\Dake\Princeton\Research\PlagborgMoller_LPVAR\MATLAB_file\Codes';
% path = '/home/dakel/Codes';
path = '/Users/ckwolf/Dropbox/Research/lp_var_simul/Codes';
cd([path '/DSGE/CKW_022020/Plots']);
addpath(genpath('../SW_G'))
addpath(genpath('../SW_MP'))

%% SETTINGS

mat_folders = {'SW_G','SW_MP'}; % Folders with files
mat_files = {{'SW_G_IV_CKW', 'SW_G_ObsShock_CKW', 'SW_G_Recursive_CKW'}, ...
             {'SW_MP_IV_CKW', 'SW_MP_ObsShock_CKW', 'SW_MP_Recursive_CKW'}}; % Files corresponding to above folders
         
output_suffix = 'png';  % File suffix

%% PLOTS

for nf=1:length(mat_folders) % For each folder...

    for ne=1:length(mat_files{nf}) % For each experiment in folder...
        
        % file/folder names
        
        exper_name = mat_files{nf}{ne}; % Name of current experiment
        file_name = strcat(exper_name); % Name of .mat results file
        output_folder = strcat(file_name, '/'); % Name of output folder
        mkdir(output_folder); % Create output folder
        
        % load results
        
        res = load(strcat(file_name, '.mat')); % Load
        settings = res.settings;
        results = res.results;
        SW_model = res.SW_model;
        clear res
        
        % compute reporting results
        
        for i_method = 1:settings.est.n_methods
    
            thisMethod = settings.est.methods_name{i_method};
            benchMethod = settings.est.methods_name{1};

            eval(['results.MSE.' thisMethod '= squeeze(mean((results.irf.' thisMethod ...
                ' - permute(SW_model.target_irf,[1 3 2])).^2, 2));']);

            eval(['results.BIAS2.' thisMethod '= (squeeze(mean(results.irf.' thisMethod ...
                ', 2)) - SW_model.target_irf).^2;']);

            eval(['results.VCE.' thisMethod '= squeeze(var(results.irf.' thisMethod ', 0, 2));']);

            eval(['results.BIAS.' thisMethod '= sqrt((squeeze(mean(results.irf.' thisMethod ...
                ', 2)) - SW_model.target_irf).^2);']);

            eval(['results.SD.' thisMethod '= sqrt(squeeze(var(results.irf.' thisMethod ', 0, 2)));']);

            eval(['results.BIASrel.' thisMethod '= sqrt((squeeze(mean(results.irf.' thisMethod ...
                ', 2)) - SW_model.target_irf).^2) ./ sqrt(mean(SW_model.target_irf.^2));']);

            eval(['results.SDrel.' thisMethod '= sqrt(squeeze(var(results.irf.' thisMethod ...
                ', 0, 2))) ./ sqrt(mean(SW_model.target_irf.^2));']);

            eval(['results.MSErel.' thisMethod '= squeeze(mean((results.irf.' thisMethod ...
                ' - permute(SW_model.target_irf,[1 3 2])).^2, 2)) ./ sqrt(mean(SW_model.target_irf.^2));']);
    
        end
        
        % plot results
        
        % bias
        
        figure
        for i_method = 1:settings.est.n_methods
            thisMethod = settings.est.methods_name{i_method};
            resultsthisMethod = eval(['results.BIASrel.' thisMethod '']);
            plot(settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Bias'), 'Interpreter', 'none');
        legend(settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'_bias'), '.', output_suffix));
%         saveas(gcf,strcat(strcat(output_folder/exper_name,'_bias'), '.', output_suffix));
        
        % std
       
        figure
        for i_method = 1:settings.est.n_methods
            thisMethod = settings.est.methods_name{i_method};
            resultsthisMethod = eval(['results.SDrel.' thisMethod '']);
            plot(settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Std'), 'Interpreter', 'none');
        legend(settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'_std'), '.', output_suffix));
%         saveas(gcf, strcat(strcat(exper_name,'_std'), '.', output_suffix));
        
        % rmse
        
        figure
        for i_method = 1:settings.est.n_methods
            thisMethod = settings.est.methods_name{i_method};
            resultsthisMethod = eval(['results.MSErel.' thisMethod '']);
            plot(settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': RMSE'), 'Interpreter', 'none');
        legend(settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'_rmse'), '.', output_suffix));
%         saveas(gcf, strcat(strcat(exper_name,'_rmse'), '.', output_suffix));
        
    end
    
end