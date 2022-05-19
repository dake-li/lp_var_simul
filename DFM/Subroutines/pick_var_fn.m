function specifications = pick_var_fn(model, settings, spec_id)
% Function for randomly choosing DGPs from the encompassing model and
% choosing IV persistence from the pre-specified grid
    % Warning: a "specification" here is a "DGP" in the paper

    % prepare
    
    specifications = settings.specifications;
    
    n_y = model.n_y;
    
    with_IV = settings.est.with_IV;
    
    if with_IV == 1
        rho = model.IV.rho;
        rho_grid = model.IV.rho_grid;
    end
    
    manual_var_select = specifications.manual_var_select;
    random_select = specifications.random_select;
    
    % randomly choose or not
    
    if random_select == 0

        % manually select

        var_select = manual_var_select;
        
        % draw IV persistency

        if with_IV == 1
            rng(spec_id, 'twister');
            rho_select_grid_idx = randi(length(rho_grid)) * ones(size(var_select, 1), 1); % arbitrarily pick one if there's multiple IV persitence levels
            rho_select = reshape(rho_grid(rho_select_grid_idx), [],1);
        end

    else

        % randomly select
        
        random_n_spec         = specifications.random_n_spec;
        random_n_var          = specifications.random_n_var;
        random_fixed_var      = specifications.random_fixed_var;
        random_fixed_pos      = specifications.random_fixed_pos;
        random_category_range = specifications.random_category_range;
        random_category_setup = specifications.random_category_setup;
        random_from_key_series = settings.specifications.random_from_key_series;
        random_exhaust_key_series = settings.specifications.random_exhaust_key_series;
        key_series = settings.specifications.random_key_series;

        if random_from_key_series == 0

            % randomly draw from the entire DFM list

            rng(spec_id, 'twister');

            var_select = nan(random_n_spec, random_n_var);
            var_select(:, 1) = random_fixed_var;
        
            for i_spec = 1:random_n_spec % draw for each specification
                
                for iv = 1:(random_n_var - 1) % draw non-fixed variables from other categories
                    
                    % check if draw from certain categories
    
                    if iv > length(random_category_setup)
                        feasible_range = 1:n_y;
                    else
                        feasible_range = [];
                        for icat = random_category_setup{iv}
                            feasible_range = [feasible_range, random_category_range(icat,1):random_category_range(icat,2)];
                        end
                    end
    
                    % draw this random variable
    
                    while true
                        
                        drawn_pos = randi(length(feasible_range));
                        drawn_ID = feasible_range(drawn_pos);
                        
                        if all(var_select(i_spec, 1:iv)~=drawn_ID) % add this new variable only if no same variable exists
                            var_select(i_spec, iv + 1) = drawn_ID;
                            break;                           
                        else % same variable exist, remove and redraw
                            feasible_range = feasible_range([1:(drawn_pos-1),(drawn_pos+1):end]);
                        end
                       
                    end
    
                end
                
                % randomly permute non-fixed variables
    
                var_select(i_spec, 2:random_n_var) = var_select(i_spec, 1 + randperm(random_n_var - 1));
                
            end

        else

            % randomly draw from some key series in the DFM list

            key_series = key_series(key_series ~= random_fixed_var);
            key_series = sort(key_series);
            key_series_category = arrayfun(@(x) sum(x > random_category_range(:,2)) + 1, key_series);

            % all combinations from key series

            all_combo_idx = nchoosek(1:length(key_series), random_n_var - 1); % index of series in each combo
            all_combo_cat = reshape(key_series_category(all_combo_idx(:)), [], random_n_var - 1); % category of series in each combo

            % keep only valid combinations
            
            if isempty(random_category_setup) % no category condition to check

                valid_combo_idx = all_combo_idx;

            else % check each category condition

                valid_label_each_cat_setup = zeros(size(all_combo_idx, 1), length(random_category_setup));
    
                for i_spec = 1:size(all_combo_idx, 1)
                    for i_setup = 1:length(random_category_setup)
                        
                        if any(ismember(all_combo_cat(i_spec, :), random_category_setup{i_setup}))
                            valid_label_each_cat_setup(i_spec, i_setup) = 1;
                        end
                    end
                end

                valid_label = all(logical(valid_label_each_cat_setup), 2);
                valid_combo_idx = all_combo_idx(valid_label, :);

            end

            % randomly permute series index and shuffle DGPs

            rng(0, 'twister'); % fix seed to pin down permutation and shuffling

            for i_spec = 1:size(valid_combo_idx, 1)
                valid_combo_idx(i_spec, :) =  valid_combo_idx(i_spec, randperm(random_n_var - 1));
            end

            valid_combo_idx = valid_combo_idx(randperm(size(valid_combo_idx, 1)), :);

            % construct exhaustive list of DGPs

            n_spec_exhaust = size(valid_combo_idx, 1);
            valid_combo_var = reshape(key_series(valid_combo_idx(:)), n_spec_exhaust, []);

            var_select = nan(n_spec_exhaust, random_n_var);
            var_select(:, 1) = random_fixed_var;
            var_select(:, 2:end) = valid_combo_var;

            % have only random_n_spec DGPs among all exhaustive DGPs

            if random_exhaust_key_series == 0

                if (spec_id-1) * random_n_spec >= n_spec_exhaust
                    error("spec_id is too large when exhausting DGPs")
                end

                % use spec_id to indicate which random_n_spec DGPs to draw

                var_select = var_select(((spec_id-1) * random_n_spec + 1):min(spec_id * random_n_spec, n_spec_exhaust), :);

            end

        end
        
        % put fixed variable at the fixed position

        var_select(:, [1 random_fixed_pos]) = var_select(:, [random_fixed_pos 1]); % put fixed var at fixed position
        
        % draw IV persistency

        if with_IV == 1
            if random_from_key_series == 1
                rng(spec_id, 'twister'); % revert back to specified seed for drawing IV persistence after exhaust from key series
            end
            rho_select_grid_idx = randi(length(rho_grid), [random_n_spec, 1]); % index of rho from rho_grid
            rho_select = reshape(rho_grid(rho_select_grid_idx), [],1); % value of rho
        end
        
    end
    
    % wrap up
    
    n_spec = size(var_select, 1);
    n_var  = size(var_select, 2);
    
    specifications.var_select = var_select;
    specifications.n_spec     = n_spec;
    specifications.n_var      = n_var;
    
    if with_IV == 1
        specifications.rho_select = rho_select;
        specifications.rho_select_grid_idx = rho_select_grid_idx;
    end

end