function IRF = compute_VARirfs_DSGE(model, settings)

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

var_select = settings.specifications.var_select;
n_spec = settings.specifications.n_spec;

recursive_shock_pos = settings.est.recursive_shock_pos;
IRF_hor = settings.est.IRF_hor;
IRF_response_var_pos = settings.est.IRF_response_var_pos;

IRF = NaN(IRF_hor, n_spec);

for i_spec = 1:n_spec

    A   = model.ABCD.A;
    B   = model.ABCD.B;
    C   = model.ABCD.C;
    C   = C(var_select(i_spec,:),:);
    D   = model.ABCD.D;
    D   = D(var_select(i_spec,:),:);
    n_s = size(A, 1);
    n_y = size(C, 1);
    
    %----------------------------------------------------------------
    % Innovation Representation
    %----------------------------------------------------------------

    % derive Sigma and K

    Sigma_states = eye(n_s);
    K_states = zeros(n_s,n_y);

    dist = 1;
    tol = 10^(-10);
    relax = 0.9;

    while dist >= tol
        Sigma_upd = (A - K_states*C)*Sigma_states*(A-K_states*C)' + B*B' + K_states*(D*D')*K_states' - B * D'*K_states' - K_states * D * B';
        K_upd = (A * Sigma_upd * C' + B * D') * (C * Sigma_upd * C' + D * D')^(-1);

        dist1 = max(max(abs(Sigma_upd - Sigma_states)));
        dist2 = max(max(abs(K_upd - K_states)));

        Sigma_states = relax * Sigma_states + (1-relax) * Sigma_upd;
        K_states = relax * K_states + (1-relax) * K_upd;

        dist = max(dist1, dist2);
    end

    %----------------------------------------------------------------
    % IRFs
    %----------------------------------------------------------------

    % recursive IRFs for y

    Sigma_u = C * Sigma_states * C' + D * D';
    G = chol(Sigma_u, 'lower');

    IRF_y = zeros(n_y, n_y, IRF_hor);    
    IRF_y(:,:,1) = G;

    for i = 1:(IRF_hor - 1)
        
        IRF_y(:,:,i+1) = C * A^(i - 1) * K_states * G;

    end
    
    % store IRF
    IRF_response = squeeze(IRF_y(IRF_response_var_pos, recursive_shock_pos, :));
    IRF_normalize = G(recursive_shock_pos, recursive_shock_pos);
    IRF(:, i_spec) = IRF_response / IRF_normalize;

end

end