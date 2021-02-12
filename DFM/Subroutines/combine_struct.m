function [DFM_estimate_combined, DF_model_combined, settings_combined, results_combined] = combine_struct(save_folder, file_prefix, spec_id_array, winsor_percent, quantiles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for id = spec_id_array

    clear DF_model DFM_estimate results settings
    load(fullfile(save_folder, strcat(file_prefix, '_', num2str(id))), 'DF_model', 'DFM_estimate', 'results', 'settings');

    %----------------------------------------------------------------
    % combine DFM_estimate
    %----------------------------------------------------------------
    if id == spec_id_array(1)
        DFM_estimate_combined = DFM_estimate;
    end
    
    %----------------------------------------------------------------
    % combine DF_model
    %----------------------------------------------------------------
    first_tier = fieldnames(DF_model);
    
    % iterate thru the first tier of the struct
    for i = 1:length(first_tier)
        
        % decide which dim to concatenate
        if any(strcmp(first_tier{i},{'R0_sq','LRV_Cov_tr_ratio','VAR_largest_root','frac_coef_for_large_lags','IV_strength'}))
            concatenate_dim = 1;
        elseif any(strcmp(first_tier{i},{'VAR_irf','IV_irf','target_irf'}))
            concatenate_dim = 2;
        else
            concatenate_dim = NaN;
        end
        
        % concatenate
        if id == 1
            DF_model_combined.(first_tier{i}) = DF_model.(first_tier{i});
        elseif ~isnan(concatenate_dim)
            DF_model_combined.(first_tier{i}) = cat(concatenate_dim, DF_model_combined.(first_tier{i}), DF_model.(first_tier{i})); 
        end
        
    end
    
    %----------------------------------------------------------------
    % combine settings
    %----------------------------------------------------------------
    first_tier = fieldnames(settings);
    
    % iterate thru the first tier of the struct
    for i = 1:length(first_tier)

        second_tier = fieldnames(settings.(first_tier{i}));
        
        % iterate thru the second tier of the struct
        for j = 1:length(second_tier)
            
            % decide which dim to concatenate
            if any(strcmp(second_tier{j},{'var_select'}))
                concatenate_dim = 1;
            else
                concatenate_dim = NaN;
            end
            
            % concatenate
            if id == 1
                settings_combined.(first_tier{i}).(second_tier{j}) = settings.(first_tier{i}).(second_tier{j});
            elseif ~isnan(concatenate_dim)
                settings_combined.(first_tier{i}).(second_tier{j}) = cat(concatenate_dim, settings_combined.(first_tier{i}).(second_tier{j}), settings.(first_tier{i}).(second_tier{j})); 
            end
            
        end

    end
    
    %----------------------------------------------------------------
    % combine results
    %----------------------------------------------------------------

    first_tier = fieldnames(results);
    
    % iterate thru the first tier of the struct
    for i = 1:length(first_tier)

        second_tier = fieldnames(results.(first_tier{i}));
        % iterate thru the second tier of the struct
        for j = 1:length(second_tier)
            
            % decide which dim to summarize (MC)
            if any(strcmp(first_tier{i},{'irf'}))
                summarize_dim = 2;
            elseif any(strcmp(first_tier{i},{'weight'}))
                summarize_dim = 3;
            elseif any(strcmp(first_tier{i},{'n_lags','largest_root','LM_stat','LM_pvalue','Granger_stat','Granger_pvalue','F_stat','F_pvalue','lambda'}))
                summarize_dim = 1;
            else
                summarize_dim = NaN;
            end
            
            % summarize (MC)
            if isnan(summarize_dim)
                results_summarized.(first_tier{i}).(second_tier{j}) = results.(first_tier{i}).(second_tier{j});
            elseif ~isnan(summarize_dim)
                this_cell_object = num2cell(results.(first_tier{i}).(second_tier{j}), summarize_dim);
                this_cell_summarized = cellfun(@(x) summ_stat(x, winsor_percent, quantiles), this_cell_object, 'UniformOutput', false);
                results_summarized.(first_tier{i}).(second_tier{j}) = cell2mat(this_cell_summarized);
            end
            
            % decide which dim to concatenate (spec)
            if any(strcmp(first_tier{i},{'irf'}))
                concatenate_dim = 3;
            elseif any(strcmp(first_tier{i},{'weight'}))
                concatenate_dim = 4;
            elseif any(strcmp(first_tier{i},{'n_lags','largest_root','LM_stat','LM_pvalue','Granger_stat','Granger_pvalue','F_stat','F_pvalue','lambda','MSE','BIAS2','VCE'}))
                concatenate_dim = 2;
            else
                concatenate_dim = NaN;
            end
            
            % concatenate (spec)
            if id == 1
                results_combined.(first_tier{i}).(second_tier{j}) = results_summarized.(first_tier{i}).(second_tier{j});
            elseif ~isnan(concatenate_dim)
                results_combined.(first_tier{i}).(second_tier{j}) = cat(concatenate_dim, results_combined.(first_tier{i}).(second_tier{j}), results_summarized.(first_tier{i}).(second_tier{j})); 
            end
           
        end

    end
    
end

% update count of specifications
settings_combined.specifications.random_n_spec = settings_combined.specifications.random_n_spec * length(spec_id_array);
settings_combined.specifications.n_spec = settings_combined.specifications.n_spec * length(spec_id_array);

% update settings
settings_combined.specifications.spec_id_array = spec_id_array;
settings_combined.simul.winsor_percent = winsor_percent;
settings_combined.simul.quantiles = quantiles;
settings_combined.simul.summ_stat_name = [{'mean','std','winsorized_mean','winsorized_std'},strcat('quant_', strsplit(num2str(quantiles)))];

end

