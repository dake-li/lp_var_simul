clear all;
addpath('./Subroutines');

% Test bias/variance formulas


%% Settings

% DGP
rho = 0.7;
sigma_2 = 2;
h = 4;
alpha = 5;

% Simulation settings
T = 1e6;
numrep = 1e3;
rng(20201223);


%% Theoretical IR and bias/var

ir = rho^h;
[bias_var, var_var, var_lp] = asy_bias_var(rho,sigma_2,alpha,h);


%% Simulate and estimate

ir_estims = nan(numrep,2);
Ahats = nan(numrep,2);
kappahats = nan(numrep,1);

randis = randi(2^32-1,numrep,1);

% for i=1:numrep
parfor i=1:numrep
   
    % Simulate
    rng(randis(i), 'twister');
    epss = randn(T,2).*[1 sigma_2];
    U = epss(:,1) + filter([1 alpha/sqrt(T)], 1, epss(:,2));
    Y = filter(1, [1 -rho], U);
    
    % Estimate
    [ir_estims(i,:),Ahats(i,:),kappahats(i)] = estim_var_lp(Y,epss(:,1),h);
    
    % Progress
    if mod(i,ceil(numrep/20))==0
        fprintf('%3d%s\n', round(100*i/numrep), '%');
    end
    
end


%% Compare bias/var

disp('Theoretical scaled bias [VAR LP]');
disp([bias_var 0]);

disp('Simulated scaled bias [VAR LP]');
disp(sqrt(T)*(mean(ir_estims)-ir));

disp('Theoretical scaled variance [VAR LP]');
disp([var_var var_lp]);

disp('Simulated scaled variance [VAR LP]');
disp(T*var(ir_estims));

