function [irf,w] = IRF_LP(Y,recurShock,respV,nlags,nhorizons)
% Auxiliary function for estimating IRFs using least-squares LP

nT = size(Y,1);

% error checking
if (nhorizons + nlags) >= nT
    error("Number of horizons too large! No obs in sample!")
end

irf = zeros(1,nhorizons + 1);

% go thru horizon 0 to horizon max
for h = 0:nhorizons
    [~,~,Bx,~,~,~,w1] = LP(Y,recurShock,respV,nlags,h); % h-step ahead LP
    irf(1,h + 1) = Bx;
    if h==0
        w=w1; % Store control data vector for later bias correction if desired
    end
end

end

