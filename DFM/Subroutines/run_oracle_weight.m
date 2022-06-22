%% COMPUTE ORACLE WEIGHTS TO AVERAGE MULTIPLE METHODS
% Dake Li, Mikkel Plagborg-Møller and Christian Wolf
% This version: 02/23/2021

clear all;
addpath(genpath(fullfile('..', '..', 'Auxiliary_Functions')))

%% SET UP DESTINATION FOLDER AND FILES

spec_id = 1; % specification choice set id array
dgp_type = 'G'; % Either 'G' or 'MP'
estimand_type = 'ObsShock'; % Either 'ObsShock', 'Recursive', or 'IV'
lag_type = 4; % No. of lags to impose in estimation, or NaN (meaning AIC)
mode_type = 1; % robustness check mode:
               % 1 (baseline), 2 (cumulative IRF), 
               % 3 (persistent DGP), 4 (persistent DGP with MN prior), 
               % 5 (small sample), 6 (salient series)

save_pre = fullfile('..', 'Results');

mode_list   = {'baseline', 'cumulative', 'persistent', 'persistent_BVAR_MN_prior' , 'small', 'salient'};
save_mode_dir = mode_list{mode_type}; % set up directory for robustness-check modes

if isnan(lag_type)
    save_suff = '_aic';
else
    save_suff = num2str(lag_type);
end
save_folder = fullfile(save_pre, save_mode_dir, strcat('lag', save_suff));

%% SET UP DETAILS FOR ORACLE WEIGHT

% select the range of methods
avg_methods_in_var_avg = 1; % Average across submodels in VAR model averaging? Otherwise, average across ultimate estimators

% select multiple estimators for averaging
% methods_to_average = {'svar','lp'};

% select multiple submodels for averaging
methods_to_average = 1:40;

% select object to minimize
minimize_object = 'MSE'; % Either 'MSE', 'BIAS2', or 'VCE'

% set up name for this oracle-weight-based estimator
save_field_name = 'var_avg_oracle';

% set up destination folder for output files
save_pre_new = [save_pre, '_oracle_weight'];

%% COMPUTE ORACLE WEIGHTS FOR SELECTED METHODS TO MINIMIZE OBJ. OF INTEREST
% (Infeasible in one sample. Compute oracle weights in population sense)

% load results
load(fullfile(save_folder, strcat('DFM_', dgp_type, '_', estimand_type, '_', num2str(spec_id))), 'DF_model', 'DFM_estimate', 'results', 'settings');

% collection of different methods
irf_estimate_collect = NaN(settings.est.IRF_hor, settings.simul.n_MC, settings.specifications.n_spec, length(methods_to_average));

for i_method = 1:length(methods_to_average)
    if avg_methods_in_var_avg == 1
        irf_estimate_collect(:,:,:,i_method) = results.submodel_irf.var_avg(i_method,:,:,:);
    else
        irf_estimate_collect(:,:,:,i_method) = results.irf.(methods_to_average{i_method});
    end
end

% permute dimension to put MC and methods first
irf_estimate_collect = permute(irf_estimate_collect, [2,4,1,3]); %MC*methods*horizon*spec

% placeholder for optimal weights and averaged estimates
method_avg_weight = NaN(length(methods_to_average), settings.est.IRF_hor, settings.specifications.n_spec); %methods*horizon*spec
method_avg_irf = NaN(settings.simul.n_MC, settings.est.IRF_hor, settings.specifications.n_spec); %MC*horizon*spec

% for each horizon*spec
for i_spec = 1:settings.specifications.n_spec
    for i_hor = 1:settings.est.IRF_hor
        
        % joint MSE, BIAS2 and VCE matrix
        joint_MSE = (irf_estimate_collect(:,:,i_hor,i_spec) - DF_model.target_irf(i_hor,i_spec))' * (irf_estimate_collect(:,:,i_hor,i_spec) - DF_model.target_irf(i_hor,i_spec)) / settings.simul.n_MC;
        joint_BIAS2 = (mean(irf_estimate_collect(:,:,i_hor,i_spec), 1) - DF_model.target_irf(i_hor,i_spec))' * (mean(irf_estimate_collect(:,:,i_hor,i_spec), 1) - DF_model.target_irf(i_hor,i_spec));
        joint_VCE = cov(irf_estimate_collect(:,:,i_hor,i_spec));
        
        % select object of interest
        switch minimize_object
            case 'MSE'
                minimize_obj_mat = joint_MSE;
            case 'BIAS2'
                minimize_obj_mat = joint_BIAS2;
            case 'VCE'
                minimize_obj_mat = joint_VCE;
        end
        
        % optimize weights to minimize object of interest
        w = quadprog(minimize_obj_mat,zeros(length(methods_to_average),1),[],[],ones(length(methods_to_average),1)',1,zeros(length(methods_to_average),1),ones(length(methods_to_average),1),[],optimoptions('quadprog','Display','off'));
        method_avg_weight(:,i_hor,i_spec) = w;
        
        % averaged estimates
        method_avg_irf(:,i_hor,i_spec) = irf_estimate_collect(:,:,i_hor,i_spec) * w;
    end
end

clear i_spec i_hor i_method irf_estimate_collect joint_MSE joint_BIAS2 joint_VCE minimize_obj w

% store results in struct
settings.est.avg_methods_in_var_avg = avg_methods_in_var_avg;
settings.est.methods_to_average = methods_to_average;
results.oracle_weight.(save_field_name) = method_avg_weight;
results.irf.(save_field_name) = permute(method_avg_irf, [2,1,3]);
if ~any(strcmp(settings.est.methods_name, save_field_name))
    settings.est.methods_name{end+1} = save_field_name;
end
settings.est.n_methods = length(settings.est.methods_name);

% re-run IRF performance summary
[results.MSE, results.BIAS2, results.VCE] = irf_perform_summary(results.irf, DF_model.target_irf, settings);

% save results
save_folder_new = fullfile(save_pre_new, strcat('lag', save_suff));
mkdir(save_folder_new);
save(fullfile(save_folder_new, strcat('DFM_', dgp_type, '_', estimand_type, '_', num2str(spec_id))), ...
    'DFM_estimate','DF_model','settings','results',...
    'spec_id','dgp_type','estimand_type','lag_type','-v7.3');

clear save_folder save_pre save_mode_dir mode_list save_suff save_folder_new save_pre_new save_field_name avg_methods_in_var_avg methods_to_average method_avg_irf method_avg_weight minimize_*