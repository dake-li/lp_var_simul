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

% colormap
n = 200;
clear cmap
cmap(1,:) = [0 0 0];
cmap(2,:) = [0.5 0.5 0.5];
cmap(3,:) = [1 1 1];

[X,Y] = meshgrid([1:3],[1:50]);

cmap = interp2(X([1,25,50],:),Y([1,25,50],:),cmap,X,Y);

% weight grid
n_weight    = 1001;
weight_grid = linspace(1,0,n_weight)';

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

            eval(['res.results.MSE.' thisMethod '= squeeze(mean((res.results.irf.' thisMethod ...
                ' - permute(res.DF_model.target_irf,[1 3 2])).^2, 2));']);

            eval(['res.results.BIAS2.' thisMethod '= (squeeze(mean(res.results.irf.' thisMethod ...
                ', 2)) - res.DF_model.target_irf).^2;']);

            eval(['res.results.VCE.' thisMethod '= squeeze(var(res.results.irf.' thisMethod ', 0, 2));']);

            eval(['res.results.BIAS.' thisMethod '= sqrt((squeeze(mean(res.results.irf.' thisMethod ...
                ', 2)) - res.DF_model.target_irf).^2);']);

            eval(['res.results.SD.' thisMethod '= sqrt(squeeze(var(res.results.irf.' thisMethod ', 0, 2)));']);

            eval(['res.results.BIASrel.' thisMethod '= sqrt((squeeze(mean(res.results.irf.' thisMethod ...
                ', 2)) - res.DF_model.target_irf).^2) ./ sqrt(mean(res.DF_model.target_irf.^2));']);

            eval(['res.results.SDrel.' thisMethod '= sqrt(squeeze(var(res.results.irf.' thisMethod ...
                ', 0, 2))) ./ sqrt(mean(res.DF_model.target_irf.^2));']);

%             eval(['res.results.MSErel.' thisMethod '= median(squeeze(mean((res.results.irf.' thisMethod ...
%                 ' - permute(res.DF_model.target_irf,[1 3 2])).^2, 2)) ./ sqrt(mean(res.DF_model.target_irf.^2)),2);']);

            eval(['res.results.MSErel.' thisMethod '= 0.5 * res.results.BIASrel.' thisMethod ...
                '.^2 + 0.5 * res.results.SDrel.' thisMethod '.^2;']);
            
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
            resultsthisMethod = median(eval(['res.results.BIASrel.' thisMethod '']),2);
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
            resultsthisMethod = median(eval(['res.results.SDrel.' thisMethod '']),2);
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
            resultsthisMethod = median(eval(['res.results.MSErel.' thisMethod '']),2);
            plot(res.settings.est.IRF_select-1,resultsthisMethod,'Linewidth',3.5)
            hold on
        end
        title(strcat(exper_name,': RMSE'), 'Interpreter', 'none');
        legend(res.settings.est.methods_name, 'Location', 'southoutside', 'Interpreter', 'none','Orientation','horizontal');
        hold off
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_rmse'), '.', output_suffix));
        
        % LP vs. VAR
        
        figure;
        
        pref_VAR = zeros(n_weight,max(res.settings.est.IRF_select),size(res.results.BIASrel.svar,2));
        for i_weight = 1:n_weight
            loss_VAR   = weight_grid(i_weight) * res.results.BIASrel.svar.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.svar.^2;
            loss_other = weight_grid(i_weight) * res.results.BIASrel.lp.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.lp.^2;

            pref_VAR(i_weight,:,:) = (loss_VAR <= loss_other);
        end
        pref_VAR = mean(pref_VAR,3);
        
        subplot(1,2,1)
        h_VAR = heatmap(pref_VAR,'ColorbarVisible','off','GridVisible','off');
        h_VAR.XLabel = 'Horizon';
        h_VAR.YLabel = 'Bias Weight';
        h_VAR.FontSize = 10;
        h_VAR.Colormap = cmap;
        h_VAR.ColorLimits = [0 1];
        h_VAR.CellLabelColor = 'none';
        CustomXLabels = string(res.settings.est.IRF_select-1);
        CustomXLabels(mod(res.settings.est.IRF_select-1,2) ~= 0) = " ";
        h_VAR.XDisplayLabels = CustomXLabels';
        CustomYLabels = weight_grid;
        CustomYLabels(mod(round(weight_grid,5),0.1) ~= 0) = " ";
        h_VAR.YDisplayLabels = CustomYLabels';
        set(struct(h_VAR).NodeChildren(3), 'XTickLabelRotation', 0);
        a2 = axes('Position', h_VAR.Position);               %new axis on top
        a2.Color = 'none';  
        a2.Title = title(strcat(exper_name,': Indifference: VAR vs. LP'), 'Interpreter', 'none');
        a2.XTick = [];
        a2.YTick = [];
        
        pref_VAR = zeros(n_weight,max(res.settings.est.IRF_select),size(res.results.BIASrel.svar,2));
        for i_weight = 1:n_weight
            loss_VAR   = weight_grid(i_weight) * res.results.BIASrel.svar.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.svar.^2;
            loss_other = weight_grid(i_weight) * res.results.BIASrel.lp_penalize.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.lp_penalize.^2;

            pref_VAR(i_weight,:,:) = (loss_VAR <= loss_other);
        end
        pref_VAR = mean(pref_VAR,3);
        
        subplot(1,2,2)
        h_VAR = heatmap(pref_VAR,'ColorbarVisible','off','GridVisible','off');
        h_VAR.XLabel = 'Horizon';
        h_VAR.YLabel = 'Bias Weight';
        h_VAR.FontSize = 10;
        h_VAR.Colormap = cmap;
        h_VAR.ColorLimits = [0 1];
        h_VAR.CellLabelColor = 'none';
        CustomXLabels = string(res.settings.est.IRF_select-1);
        CustomXLabels(mod(res.settings.est.IRF_select-1,2) ~= 0) = " ";
        h_VAR.XDisplayLabels = CustomXLabels';
        CustomYLabels = weight_grid;
        CustomYLabels(mod(round(weight_grid,5),0.1) ~= 0) = " ";
        h_VAR.YDisplayLabels = CustomYLabels';
        set(struct(h_VAR).NodeChildren(3), 'XTickLabelRotation', 0);
        a2 = axes('Position', h_VAR.Position);               %new axis on top
        a2.Color = 'none';  
        a2.Title = title(strcat(exper_name,': Indifference: VAR vs. pen. LP'), 'Interpreter', 'none');
        a2.XTick = [];
        a2.YTick = [];
        
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1) pos(2) 1.6*pos(3) 0.8*pos(4)]);
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_varvslp'), '.', output_suffix));
        
        % VAR methods
        
        figure;
        
        pref_VAR = zeros(n_weight,max(res.settings.est.IRF_select),size(res.results.BIASrel.svar,2));
        for i_weight = 1:n_weight
            loss_VAR   = weight_grid(i_weight) * res.results.BIASrel.svar.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.svar.^2;
            loss_other = weight_grid(i_weight) * res.results.BIASrel.svar_corrbias.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.svar_corrbias.^2;

            pref_VAR(i_weight,:,:) = (loss_VAR <= loss_other);
        end
        pref_VAR = mean(pref_VAR,3);
        
        subplot(1,3,1)
        h_VAR = heatmap(pref_VAR,'ColorbarVisible','off','GridVisible','off');
        h_VAR.XLabel = 'Horizon';
        h_VAR.YLabel = 'Bias Weight';
        h_VAR.FontSize = 10;
        h_VAR.Colormap = cmap;
        h_VAR.ColorLimits = [0 1];
        h_VAR.CellLabelColor = 'none';
        CustomXLabels = string(res.settings.est.IRF_select-1);
        CustomXLabels(mod(res.settings.est.IRF_select-1,2) ~= 0) = " ";
        h_VAR.XDisplayLabels = CustomXLabels';
        CustomYLabels = weight_grid;
        CustomYLabels(mod(round(weight_grid,5),0.1) ~= 0) = " ";
        h_VAR.YDisplayLabels = CustomYLabels';
        set(struct(h_VAR).NodeChildren(3), 'XTickLabelRotation', 0);
        a2 = axes('Position', h_VAR.Position);               %new axis on top
        a2.Color = 'none';  
        a2.Title = title(strcat(exper_name,': Indifference: VAR vs. bias-c. VAR'), 'Interpreter', 'none');
        a2.XTick = [];
        a2.YTick = [];
        
        pref_VAR = zeros(n_weight,max(res.settings.est.IRF_select),size(res.results.BIASrel.svar,2));
        for i_weight = 1:n_weight
            loss_VAR   = weight_grid(i_weight) * res.results.BIASrel.svar.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.svar.^2;
            loss_other = weight_grid(i_weight) * res.results.BIASrel.bvar.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.bvar.^2;

            pref_VAR(i_weight,:,:) = (loss_VAR <= loss_other);
        end
        pref_VAR = mean(pref_VAR,3);
        
        subplot(1,3,2)
        h_VAR = heatmap(pref_VAR,'ColorbarVisible','off','GridVisible','off');
        h_VAR.XLabel = 'Horizon';
        h_VAR.YLabel = 'Bias Weight';
        h_VAR.FontSize = 10;
        h_VAR.Colormap = cmap;
        h_VAR.ColorLimits = [0 1];
        h_VAR.CellLabelColor = 'none';
        CustomXLabels = string(res.settings.est.IRF_select-1);
        CustomXLabels(mod(res.settings.est.IRF_select-1,2) ~= 0) = " ";
        h_VAR.XDisplayLabels = CustomXLabels';
        CustomYLabels = weight_grid;
        CustomYLabels(mod(round(weight_grid,5),0.1) ~= 0) = " ";
        h_VAR.YDisplayLabels = CustomYLabels';
        set(struct(h_VAR).NodeChildren(3), 'XTickLabelRotation', 0);
        a2 = axes('Position', h_VAR.Position);               %new axis on top
        a2.Color = 'none';  
        a2.Title = title(strcat(exper_name,': Indifference: VAR vs. BVAR'), 'Interpreter', 'none');
        a2.XTick = [];
        a2.YTick = [];
        
        pref_VAR = zeros(n_weight,max(res.settings.est.IRF_select),size(res.results.BIASrel.svar,2));
        for i_weight = 1:n_weight
            loss_VAR   = weight_grid(i_weight) * res.results.BIASrel.svar.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.svar.^2;
            loss_other = weight_grid(i_weight) * res.results.BIASrel.var_avg.^2 + (1-weight_grid(i_weight)) * res.results.SDrel.var_avg.^2;

            pref_VAR(i_weight,:,:) = (loss_VAR <= loss_other);
        end
        pref_VAR = mean(pref_VAR,3);
        
        subplot(1,3,3)
        h_VAR = heatmap(pref_VAR,'ColorbarVisible','off','GridVisible','off');
        h_VAR.XLabel = 'Horizon';
        h_VAR.YLabel = 'Bias Weight';
        h_VAR.FontSize = 10;
        h_VAR.Colormap = cmap;
        h_VAR.ColorLimits = [0 1];
        h_VAR.CellLabelColor = 'none';
        CustomXLabels = string(res.settings.est.IRF_select-1);
        CustomXLabels(mod(res.settings.est.IRF_select-1,2) ~= 0) = " ";
        h_VAR.XDisplayLabels = CustomXLabels';
        CustomYLabels = weight_grid;
        CustomYLabels(mod(round(weight_grid,5),0.1) ~= 0) = " ";
        h_VAR.YDisplayLabels = CustomYLabels';
        set(struct(h_VAR).NodeChildren(3), 'XTickLabelRotation', 0);
        a2 = axes('Position', h_VAR.Position);               %new axis on top
        a2.Color = 'none';  
        a2.Title = title(strcat(exper_name,': Indifference: VAR vs. avg. VAR'), 'Interpreter', 'none');
        a2.XTick = [];
        a2.YTick = [];
        
        pos = get(gcf, 'Position');
        set(gcf, 'Position', [pos(1) pos(2) 2*pos(3) 0.8*pos(4)]);
        set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf, strcat(strcat(output_folder,'/',exper_name,'_vars_comp'), '.', output_suffix));
        
    end
    
end