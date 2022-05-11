function specifications = pick_var_fn(model, settings, spec_id)
% Function for randomly choosing DGPs from the encompassing model and
% choosing IV persistence from the pre-specified grid
    % Warning: a "specification" here is a "DGP" in the paper

    % prepare
    
    rng(spec_id, 'twister');

    specifications = settings.specifications;
    
    n_y = model.n_y;
    
    with_IV = settings.est.with_IV;
    
    if with_IV == 1
        rho = model.IV.rho;
        rho_grid = model.IV.rho_grid;
    end
    
    manual_var_select = specifications.manual_var_select;
    random_select = specifications.random_select;
    random_category_range = specifications.random_category_range;

    % randomly choose
    
    if random_select == 0

        var_select = manual_var_select;
        
        % draw IV persistency
        if with_IV == 1
            rho_select_grid_idx = randi(length(rho_grid)) * ones(size(var_select, 1), 1); % arbitrarily pick one if there's multiple IV persitence levels
            rho_select = reshape(rho_grid(rho_select_grid_idx), [],1);
        end

    else
        
        random_n_spec         = specifications.random_n_spec;
        random_n_var          = specifications.random_n_var;
        random_fixed_var      = specifications.random_fixed_var;
        random_fixed_pos      = specifications.random_fixed_pos;
        random_category_setup = specifications.random_category_setup;
        
        var_select = nan(random_n_spec, random_n_var);
        var_select(:, 1) = random_fixed_var;
        fixed_category = sum(random_fixed_var > random_category_range(:,2)) + 1; % category of fixed variable
        
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
                    end
                   
                end

            end
            
            % randomly permute non-fixed variables
            var_select(i_spec, 2:random_n_var) = var_select(i_spec, 1 + randperm(random_n_var - 1));
            
        end
        
        % put fixed variable at the fixed position
        var_select(:, [1 random_fixed_pos]) = var_select(:, [random_fixed_pos 1]); % put fixed var at fixed position
        
        % draw IV persistency
        if with_IV == 1
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