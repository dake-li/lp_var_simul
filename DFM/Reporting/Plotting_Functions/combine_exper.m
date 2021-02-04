function res = combine_exper(res, res_part)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%----------------------------------------------------------------
% combine DFM_estimate
%----------------------------------------------------------------
% res.DFM_estimate = res.DFM_estimate;

%----------------------------------------------------------------
% combine DF_model
%----------------------------------------------------------------
first_tier = fieldnames(res.DF_model);

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
    if ~isnan(concatenate_dim)
        res.DF_model.(first_tier{i}) = cat(concatenate_dim, res.DF_model.(first_tier{i}), res_part.DF_model.(first_tier{i})); 
    end

end

%----------------------------------------------------------------
% combine settings
%----------------------------------------------------------------
first_tier = fieldnames(res.settings);

% iterate thru the first tier of the struct
for i = 1:length(first_tier)

    second_tier = fieldnames(res.settings.(first_tier{i}));

    % iterate thru the second tier of the struct
    for j = 1:length(second_tier)

        % decide which dim to concatenate
        if any(strcmp(second_tier{j},{'var_select'}))
            concatenate_dim = 1;
        else
            concatenate_dim = NaN;
        end

        % concatenate
        if ~isnan(concatenate_dim)
            res.settings.(first_tier{i}).(second_tier{j}) = cat(concatenate_dim, res.settings.(first_tier{i}).(second_tier{j}), res_part.settings.(first_tier{i}).(second_tier{j})); 
        end

    end

end

%----------------------------------------------------------------
% combine results
%----------------------------------------------------------------

first_tier = fieldnames(res.results);

% iterate thru the first tier of the struct
for i = 1:length(first_tier)

    second_tier = fieldnames(res.results.(first_tier{i}));
    % iterate thru the second tier of the struct
    for j = 1:length(second_tier)

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
        if ~isnan(concatenate_dim)
            res.results.(first_tier{i}).(second_tier{j}) = cat(concatenate_dim, res.results.(first_tier{i}).(second_tier{j}), res_part.results.(first_tier{i}).(second_tier{j})); 
        end

    end

end

% update count of specifications
res.settings.specifications.random_n_spec = size(res.DF_model.target_irf, 2);
res.settings.specifications.n_spec = size(res.DF_model.target_irf, 2);

end

