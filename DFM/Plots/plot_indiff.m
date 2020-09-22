%% PLOT SIMULATION RESULTS
% this version: 30/8/2020

%% HOUSEKEEPING

clc
clear all
close all

path = '/Users/christianwolf/Documents/GitHub/lp_var_simul/Reporting';
cd(path);

%% SETTINGS

mat_rootfolder = '/Users/christianwolf/Dropbox/Research/lp_var_simul/stored_result/small_run_0916'; % Root folder with files
% mat_folders = {'lag4/Results', 'lag8/Results'}; % Folders with files
mat_folders = {'lag4/Results'}; % Folders with files
% mat_files = {'DFM_G_IV', 'DFM_G_ObsShock', 'DFM_G_Recursive', ...
%              'DFM_MP_IV', 'DFM_MP_ObsShock', 'DFM_MP_Recursive'}; % Files in each of the above folders
% mat_files = {'DFM_MP_Recursive'}; % Files in each of the above folders
mat_files = {'DFM_G_ObsShock', 'DFM_G_Recursive', ...
             'DFM_MP_ObsShock', 'DFM_MP_Recursive'}; % Files in each of the above folders

% Output settings for figures
output_dir = 'fig';     % Folder
output_suffix = 'png';  % File suffix

%% RAW FIGURES

for nf=1:length(mat_folders) % For each folder...

    for ne=1:length(mat_files) % For each experiment in folder...
        
        % file/folder names
        
        exper_name = mat_files{ne}; % Name of current experiment
        file_name = fullfile(mat_folders{nf}, exper_name); % Name of .mat results file
        output_folder = fullfile(output_dir, file_name); % Name of output folder
        mkdir(output_folder); % Create output folder        

        % load results

        res = load(fullfile(mat_rootfolder, strcat(file_name, '.mat'))); % Load
        horzs = res.settings.est.IRF_select; % Impulse response horizons
        methods_name = res.settings.est.methods_name; % Names of estimation methods

        % compute reporting results
        
        for i_method = 1:res.settings.est.n_methods
    
            thisMethod = res.settings.est.methods_name{i_method};
            benchMethod = res.settings.est.methods_name{1};

            eval(['res.results.MSE.' thisMethod '= median(squeeze(mean((res.results.irf.' thisMethod ...
                ' - permute(res.DF_model.target_irf,[1 3 2])).^2, 2)),2);']);

            eval(['res.results.BIAS2.' thisMethod '= median((squeeze(mean(res.results.irf.' thisMethod ...
                ', 2)) - res.DF_model.target_irf).^2,2);']);

            eval(['res.results.VCE.' thisMethod '= median(squeeze(var(res.results.irf.' thisMethod ', 0, 2)),2);']);

            eval(['res.results.BIAS.' thisMethod '= median(sqrt((squeeze(mean(res.results.irf.' thisMethod ...
                ', 2)) - res.DF_model.target_irf).^2),2);']);

            eval(['res.results.SD.' thisMethod '= median(sqrt(squeeze(var(res.results.irf.' thisMethod ', 0, 2))),2);']);

            eval(['res.results.BIASrel.' thisMethod '= median(sqrt((squeeze(mean(res.results.irf.' thisMethod ...
                ', 2)) - res.DF_model.target_irf).^2) ./ sqrt(mean(res.DF_model.target_irf.^2)),2);']);

            eval(['res.results.SDrel.' thisMethod '= median(sqrt(squeeze(var(res.results.irf.' thisMethod ...
                ', 0, 2))) ./ sqrt(mean(res.DF_model.target_irf.^2)),2);']);

            eval(['res.results.MSErel.' thisMethod '= median(squeeze(mean((res.results.irf.' thisMethod ...
                ' - permute(res.DF_model.target_irf,[1 3 2])).^2, 2)) ./ sqrt(mean(res.DF_model.target_irf.^2)),2);']);
            
            eval(['res.results.biasweight.' thisMethod ...
                '= (res.results.SDrel.' thisMethod '.^2 - res.results.SDrel.' benchMethod ...
                '.^2 ) ./ (res.results.BIASrel.' benchMethod '.^2 - res.results.BIASrel.' thisMethod ...
                '.^2 + res.results.SDrel.' thisMethod '.^2 - res.results.SDrel.' benchMethod '.^2);']);
            eval(['res.results.biasweight.' thisMethod '(res.results.biasweight.' thisMethod '>1) = NaN;']);
            eval(['res.results.biasweight.' thisMethod '(res.results.biasweight.' thisMethod '<0) = NaN;']);
    
        end
        
        % plot results
        
        % bias
        
        figure;
        for i_method = 1:res.settings.est.n_methods
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.BIASrel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Bias'), 'Interpreter', 'none');
        legend(res.settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_bias'), '.', output_suffix));
        
        % std
       
        figure;
        for i_method = 1:res.settings.est.n_methods
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.SDrel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Std'), 'Interpreter', 'none');
        legend(res.settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_std'), '.', output_suffix));
        
        % rmse
        
        figure;
        for i_method = 1:res.settings.est.n_methods
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.MSErel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': RMSE'), 'Interpreter', 'none');
        legend(res.settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_rmse'), '.', output_suffix));
        
        % LP vs. VAR
        
        figure;
        subplot(2,3,1)
        for i_method = [1,4,5]
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.BIASrel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Bias'), 'Interpreter', 'none');
        hold off
        
        subplot(2,3,2)
        for i_method = [1,4,5]
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.SDrel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Std'), 'Interpreter', 'none');
        legend(res.settings.est.methods_name{[1,4,5]}, 'Location', 'north', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        
        subplot(2,3,3)
        for i_method = [1,4,5]
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.MSErel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': RMSE'), 'Interpreter', 'none');
        hold off
        
        pref_VAR = res.results.BIASrel.svar < res.results.BIASrel.lp & res.results.SDrel.svar < res.results.SDrel.lp;
        pref_LP  = res.results.BIASrel.svar > res.results.BIASrel.lp & res.results.SDrel.svar > res.results.SDrel.lp;
        
        subplot(2,3,4)
        for i = 1:max(res.settings.est.IRF_select-1)
            if pref_VAR(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[135/255 206/255 250/255],'EdgeColor','none')
            elseif pref_LP(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[255/255 204/255 203/255],'EdgeColor','none')
            end
            hold on
        end
        hold on
        thisMethod = res.settings.est.methods_name{4};
        resultsthisMethod = eval(['res.results.biasweight.' thisMethod '']);
        plot(res.settings.est.IRF_select-1,resultsthisMethod,'-s',...
            'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'Linewidth',3,'Color',[0 0 0])
        hold on
        title(strcat(exper_name,': Indifference: VAR vs. LP'), 'Interpreter', 'none');
        hold off
        
        pref_VAR = res.results.BIASrel.svar < res.results.BIASrel.lp_penalize & res.results.SDrel.svar < res.results.SDrel.lp_penalize;
        pref_LP  = res.results.BIASrel.svar > res.results.BIASrel.lp_penalize & res.results.SDrel.svar > res.results.SDrel.lp_penalize;
        
        subplot(2,3,5)
        for i = 1:max(res.settings.est.IRF_select)-1
            if pref_VAR(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[135/255 206/255 250/255],'EdgeColor','none')
            elseif pref_LP(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[255/255 204/255 203/255],'EdgeColor','none')
            end
            hold on
        end
        thisMethod = res.settings.est.methods_name{5};
        resultsthisMethod = eval(['res.results.biasweight.' thisMethod '']);
        plot(res.settings.est.IRF_select-1,resultsthisMethod,'-s',...
            'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'Linewidth',3,'Color',[0 0 0])
        hold on
        title(strcat(exper_name,': Indifference: VAR vs. pen. LP'), 'Interpreter', 'none');
        hold off
        
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1) pos(2) 1.75*pos(3) 1.5*pos(4)]);
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_varvslp'), '.', output_suffix));
        
        % VAR methods
        
        figure;
        subplot(2,3,1)
        for i_method = [1,2,3]
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.BIASrel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Bias'), 'Interpreter', 'none');
        hold off
        
        subplot(2,3,2)
        for i_method = [1,2,3]
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.SDrel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': Std'), 'Interpreter', 'none');
        legend(res.settings.est.methods_name{[1,2,3,6]}, 'Location', 'north', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        
        subplot(2,3,3)
        for i_method = [1,2,3]
            thisMethod = res.settings.est.methods_name{i_method};
            resultsthisMethod = eval(['res.results.MSErel.' thisMethod '']);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': RMSE'), 'Interpreter', 'none');
        hold off
        
        pref_VAR   = res.results.BIASrel.svar < res.results.BIASrel.svar_corrbias & res.results.SDrel.svar < res.results.SDrel.svar_corrbias;
        pref_cbias = res.results.BIASrel.svar > res.results.BIASrel.svar_corrbias & res.results.SDrel.svar > res.results.SDrel.svar_corrbias;
        
        subplot(2,3,4)
        for i = 1:max(res.settings.est.IRF_select-1)
            if pref_VAR(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[135/255 206/255 250/255],'EdgeColor','none')
            elseif pref_cbias(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[255/255 204/255 203/255],'EdgeColor','none')
            end
            hold on
        end
        hold on
        thisMethod = res.settings.est.methods_name{2};
        resultsthisMethod = eval(['res.results.biasweight.' thisMethod '']);
        plot(res.settings.est.IRF_select-1,resultsthisMethod,'-s',...
            'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'Linewidth',3,'Color',[0 0 0])
        hold on
        title(strcat(exper_name,': Indifference: VAR vs. bias-c. VAR'), 'Interpreter', 'none');
        hold off
        
        pref_VAR   = res.results.BIASrel.svar < res.results.BIASrel.bvar & res.results.SDrel.svar < res.results.SDrel.bvar;
        pref_bvar  = res.results.BIASrel.svar > res.results.BIASrel.bvar & res.results.SDrel.svar > res.results.SDrel.bvar;
        
        subplot(2,3,5)
        for i = 1:max(res.settings.est.IRF_select-1)
            if pref_VAR(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[135/255 206/255 250/255],'EdgeColor','none')
            elseif pref_bvar(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[255/255 204/255 203/255],'EdgeColor','none')
            end
            hold on
        end
        hold on  
        thisMethod = res.settings.est.methods_name{3};
        resultsthisMethod = eval(['res.results.biasweight.' thisMethod '']);
        plot(res.settings.est.IRF_select-1,resultsthisMethod,'-s',...
            'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'Linewidth',3,'Color',[0 0 0])
        hold on
        title(strcat(exper_name,': Indifference: VAR vs. BVAR'), 'Interpreter', 'none');
        hold off
        
        pref_VAR    = res.results.BIASrel.svar < res.results.BIASrel.var_avg & res.results.SDrel.svar < res.results.SDrel.var_avg;
        pref_VARavg = res.results.BIASrel.svar > res.results.BIASrel.var_avg & res.results.SDrel.svar > res.results.SDrel.var_avg;
        
        subplot(2,3,6)
        for i = 1:max(res.settings.est.IRF_select-1)
            if pref_VAR(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[135/255 206/255 250/255],'EdgeColor','none')
            elseif pref_VARavg(i) == 1
                patch([i-1 i i i-1],[0 0 1 1],[255/255 204/255 203/255],'EdgeColor','none')
            end
            hold on
        end
        hold on   
        thisMethod = res.settings.est.methods_name{6};
        resultsthisMethod = eval(['res.results.biasweight.' thisMethod '']);
        plot(res.settings.est.IRF_select-1,resultsthisMethod,'-s',...
            'MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'Linewidth',3,'Color',[0 0 0])
        hold on
        title(strcat(exper_name,': Indifference: VAR vs. avg. VAR'), 'Interpreter', 'none');
        hold off
        
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1) pos(2) 1.75*pos(3) 1.5*pos(4)]);
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_vars_comp'), '.', output_suffix));
        
    end
    
end