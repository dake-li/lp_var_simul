function x_cumul = cumsum_nan(x)

    [n,m] = size(x);

    % Cumulate x along rows, accounting for NaNs
    first_nm = find(~isnan(x(:,1)),1); % First non-missing obs. in first column
    x_cumul = nan(n+1,m);
    x_cumul(first_nm:end,:) = [zeros(1,m);cumsum(x(first_nm:end,:),1)];

end