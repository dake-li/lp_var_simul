%% SPECIFIC SETTINGS FOR PRE-SPECIFIED ROBUSTNESS-CHECK MODES

% set up directory for robustness-check modes

mode_list   = {'baseline', 'small', 'large', 'salient', 'more', 'diff'};
save_mode_dir = mode_list{mode_type};

% rewrite some baseline settings in "shared.m" for different robustness check modes

switch mode_type

    case 1 % baseline (either levels or first differences)
        
        % rewrite nothing and use all the settings in "shared.m"

    case 2 % small sample

        % sample settings

        settings.simul.T      = 100; % time periods for each simulation

        % lag specification (accordingly reduce size of the largest model in VAR model averaging)

        settings.est.n_lags_max     = 12; % maximal lag length for info criteria

    case 3 % large sample

        % sample settings

        settings.simul.T      = 720; % time periods for each simulation

    case 4 % salient series

        % selection of DGPs from encompassing model

        settings.specifications.random_from_key_series = 1; % randomly select from some key series in DFM list?

    case 5 % more observables
        
        % selection of DGPs from encompassing model

        settings.specifications.random_n_var = 7; % larger number of variables in each specification

    case 6 % data in first differences

        DF_model.levels = 0;
        DF_model.n_lags_fac = 2; % 2 lags in factor process, as in Stock & Watson
        DF_model.n_lags_uar = 2; % 2 lags in idiosyncratic disturbance processes, as in S&W
        settings.est.prior.towards_random_walk = 0; % shrink towards white noise

end