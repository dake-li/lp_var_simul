function res = combine_struct(res, res_part, summ_option, DGP_selected)
% Function to combine the results across different DGPs
    % This function will first summarize res_part across MCs if both
    % res_part and summ_option are non-empty
    
    % It will then concatenate res and res_part along the dim of DGPs if
    % both res and res_part are non-empty
    
    % It will finally select a subset of DGPs in res based on the selection
    % label DGP_selected

%----------------------------------------------------------------
% check mode
%----------------------------------------------------------------

% concatenate different subsets of DGPs or not
if isempty(res) && isempty(res_part)
    error('No input!');
elseif isempty(res)
    mode_concat = 0;
    empty_res = 1; % if input of res is empty
elseif isempty(res_part)
    mode_concat = 0;
    empty_res = 0;
else
    mode_concat = 1;
    empty_res = 0;
end

% summarize across MCs or not
if isempty(res_part)
    mode_summ = 0;
elseif isempty(summ_option)
    mode_summ = 0;
else
    winsor_percent = summ_option.winsor_percent;
    quantiles = summ_option.quantiles;
    summ_stat_name = summ_option.summ_stat_name;
    mode_summ = 1;
end

% select a subset of DGPs or not
if isempty(DGP_selected)
    mode_select = 0;
elseif all(DGP_selected == 1)
    mode_select = 0;
elseif all(DGP_selected == 0)
    error('No DGP selected!');
else
    DGP_selected_index = find(DGP_selected == 1);
    mode_select = 1;
end

%----------------------------------------------------------------
% struct: DFM_estimate
%----------------------------------------------------------------

if empty_res == 1
    res.DFM_estimate = res_part.DFM_estimate;
end

%----------------------------------------------------------------
% struct: DF_model
%----------------------------------------------------------------

if empty_res == 1
    first_tier = fieldnames(res_part.DF_model);
else
    first_tier = fieldnames(res.DF_model);
end

% iterate thru the first tier of the struct
for i = 1:length(first_tier)

    % decide which dim to concatenate
    if any(strcmp(first_tier{i},{'R0_sq','LRV_Cov_tr_ratio','VAR_largest_root','frac_coef_for_large_lags','IV_strength'}))
        concatenate_dim = 1;
    elseif any(strcmp(first_tier{i},{'VAR_irf','normalized_irf','target_irf'}))
        concatenate_dim = 2;
    else
        concatenate_dim = NaN;
    end

    % concatenate
    if mode_concat == 0
        if empty_res == 1
            res.DF_model.(first_tier{i}) = res_part.DF_model.(first_tier{i});
        end
    else
        if ~isnan(concatenate_dim)
            res.DF_model.(first_tier{i}) = cat(concatenate_dim, res.DF_model.(first_tier{i}), res_part.DF_model.(first_tier{i})); 
        end
    end
    
    % select
    if mode_select == 1
        if ~isnan(concatenate_dim)
            all_index = arrayfun(@(x) 1:x, size(res.DF_model.(first_tier{i})), 'UniformOutput', false);
            all_index{concatenate_dim} = DGP_selected_index;
            res.DF_model.(first_tier{i}) = res.DF_model.(first_tier{i})(all_index{:});
        end
    end

end

%----------------------------------------------------------------
% struct: settings
%----------------------------------------------------------------

if empty_res == 1
    first_tier = fieldnames(res_part.settings);
else
    first_tier = fieldnames(res.settings);
end

% iterate thru the first tier of the struct
for i = 1:length(first_tier)

    if empty_res == 1
        second_tier = fieldnames(res_part.settings.(first_tier{i}));
    else
        second_tier = fieldnames(res.settings.(first_tier{i}));
    end

    % iterate thru the second tier of the struct
    for j = 1:length(second_tier)

        % decide which dim to concatenate
        if any(strcmp(second_tier{j},{'var_select','rho_select','rho_select_grid_idx'}))
            concatenate_dim = 1;
        else
            concatenate_dim = NaN;
        end

        
        % concatenate
        if mode_concat == 0
            if empty_res == 1
                res.settings.(first_tier{i}).(second_tier{j}) = res_part.settings.(first_tier{i}).(second_tier{j});
            end
        else
            if ~isnan(concatenate_dim)
                res.settings.(first_tier{i}).(second_tier{j}) = cat(concatenate_dim, res.settings.(first_tier{i}).(second_tier{j}), res_part.settings.(first_tier{i}).(second_tier{j})); 
            end
        end
        
        % select
        if mode_select == 1
            if ~isnan(concatenate_dim)
                all_index = arrayfun(@(x) 1:x, size(res.settings.(first_tier{i}).(second_tier{j})), 'UniformOutput', false);
                all_index{concatenate_dim} = DGP_selected_index;
                res.settings.(first_tier{i}).(second_tier{j}) = res.settings.(first_tier{i}).(second_tier{j})(all_index{:});
            end
        end

    end

end

%----------------------------------------------------------------
% struct: results
%----------------------------------------------------------------

if empty_res == 1
    first_tier = fieldnames(res_part.results);
else
    first_tier = fieldnames(res.results);
end

% iterate thru the first tier of the struct
for i = 1:length(first_tier)

    if empty_res == 1
        second_tier = fieldnames(res_part.results.(first_tier{i}));
    else
        second_tier = fieldnames(res.results.(first_tier{i}));
    end
    
    % iterate thru the second tier of the struct
    for j = 1:length(second_tier)

        % decide which dim to summarize (MC)
        if any(strcmp(first_tier{i},{'irf'}))
            summarize_dim = 2;
        elseif any(strcmp(first_tier{i},{'weight','submodel_irf'}))
            summarize_dim = 3;
        elseif any(strcmp(first_tier{i},{'n_lags','largest_root','LM_stat','LM_pvalue','Hausman_stat','Hausman_pvalue','Granger_stat','Granger_pvalue','F_stat','F_pvalue','lambda'}))
            summarize_dim = 1;
        else
            summarize_dim = NaN;
        end

        % summarize (MC)
        if mode_summ == 1
            if ~isnan(summarize_dim)
                this_cell_object = num2cell(res_part.results.(first_tier{i}).(second_tier{j}), summarize_dim);
                this_cell_summarized = cellfun(@(x) summ_stat(x, winsor_percent, quantiles), this_cell_object, 'UniformOutput', false);
                res_part.results.(first_tier{i}).(second_tier{j}) = cell2mat(this_cell_summarized);
            end
        end
        
        % decide which dim to concatenate (spec)
        if any(strcmp(first_tier{i},{'irf','oracle_weight'}))
            concatenate_dim = 3;
        elseif any(strcmp(first_tier{i},{'weight','submodel_irf'}))
            concatenate_dim = 4;
        elseif any(strcmp(first_tier{i},{'n_lags','largest_root','LM_stat','LM_pvalue','Hausman_stat','Hausman_pvalue','Granger_stat','Granger_pvalue','F_stat','F_pvalue','lambda','MSE','BIAS2','VCE'}))
            concatenate_dim = 2;
        else
            concatenate_dim = NaN;
        end

        % concatenate
        if mode_concat == 0
            if empty_res == 1
                res.results.(first_tier{i}).(second_tier{j}) = res_part.results.(first_tier{i}).(second_tier{j});
            end
        else
            if ~isnan(concatenate_dim)
                res.results.(first_tier{i}).(second_tier{j}) = cat(concatenate_dim, res.results.(first_tier{i}).(second_tier{j}), res_part.results.(first_tier{i}).(second_tier{j})); 
            end
        end
        
        % select
        if mode_select == 1
            if ~isnan(concatenate_dim)
                all_index = arrayfun(@(x) 1:x, size(res.results.(first_tier{i}).(second_tier{j})), 'UniformOutput', false);
                all_index{concatenate_dim} = DGP_selected_index;
                res.results.(first_tier{i}).(second_tier{j}) = res.results.(first_tier{i}).(second_tier{j})(all_index{:});
            end
        end

    end

end

%----------------------------------------------------------------
% update info
%----------------------------------------------------------------

% update count of specifications
res.settings.specifications.random_n_spec = size(res.settings.specifications.var_select,1);
res.settings.specifications.n_spec = size(res.settings.specifications.var_select,1);

% update settings of summary across MCs
if mode_summ == 1
    res.settings.simul.winsor_percent = winsor_percent;
    res.settings.simul.quantiles = quantiles;
    res.settings.simul.summ_stat_name = summ_stat_name;
end

end

