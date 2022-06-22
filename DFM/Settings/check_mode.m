%% SPECIFIC SETTINGS FOR PRE-SPECIFIED ROBUSTNESS-CHECK MODES

% set up directory for robustness-check modes

mode_list   = {'baseline', 'cumulative', 'persistent', 'persistent_BVAR_MN_prior' , 'small', 'salient'};
save_mode_dir = mode_list{mode_type};

% rewrite some baseline settings in "shared.m" for different robustness check modes

switch mode_type

    case 1 % baseline
        
        % rewrite nothing and use all the settings in "shared.m"

    case 2 % cumulative IRF
        
        % rewrite nothing and use all the settings in "shared.m"
        % cumulative IRF will be imputed in "run_combine.m"

    case 3 % persistent DGP
        
        % scale up calibrated persistency of factors

        DF_model.fac_persist.scale = 1; % scale calibrated persistence of factors?

    case 4 % persistent DGP with MN prior

        % scale up calibrated persistency of factors

        DF_model.fac_persist.scale = 1; % scale calibrated persistence of factors?

        % BVAR prior

        settings.est.prior.towards_random_walk = 1; % prior shrinking towards random walk? otherwise towards zero

    case 5 % small sample

        % sample settings

        settings.simul.T      = 100; % time periods for each simulation

        % lag specification (accordingly reduce size of the largest model in VAR model averaging)

        settings.est.n_lags_max     = 12; % maximal lag length for info criteria

    case 6 % salient series

        % selection of DGPs from encompassing model

        settings.specifications.random_from_key_series = 1; % randomly select from some key series in DFM list?

end