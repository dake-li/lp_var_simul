rho = -1 + 2 * rand;
sigma = rand;

n_h = 100;
LHS1 = NaN(n_h,1);
RHS1 = NaN(n_h,1);

for h = 1:n_h
    LHS1(h) = rho^(2*h) * (sigma^2 + 1) * (sigma^2 + rho^2) * (1-rho^2) * (1 + rho^(-2)) ...
        + rho^(2*(h-1)) * (sigma^2 + 1) * (h-1)^2 * (1-rho^2)^2;
    RHS1(h) = (sigma^2 + 1) * (sigma^2 + rho^2) * (1 - rho^(2*(h+1)));
end

min(LHS1(2:end) < RHS1(2:end))

LHS2 = NaN(n_h,1);
RHS2 = NaN(n_h,1);

for h = 1:n_h
    LHS2(h) = sigma^4 * rho^(2*(h-1)) + sigma^2 * rho^(2*h) + sigma^2 * rho^(2*(h-1)) + rho^(2*h) - sigma^4 * rho^(2*(h+1)) - sigma^2 * rho^(2*(h+2)) - sigma^2 * rho^(2*(h+1)) - rho^(2*(h+2)) ...
        + (h-1)^2 * (sigma^2 * rho^(2*(h-1)) - 2 * sigma^2 * rho^(2*h) + rho^(2*(h+1)) * sigma^2 + rho^(2*(h-1)) - 2 * rho^(2*h) + rho^(2*(h+1)));
    RHS2(h) = sigma^4 + sigma^2 * rho^2 + sigma^2 + rho^2 - sigma^4 * rho^(2*(h+1)) - sigma^2 * rho^(2*(h+2)) - sigma^2 * rho^(2*(h+1)) - rho^(2*(h+2));
end

min(LHS2(2:end) < RHS2(2:end))

LHS3 = NaN(n_h,1);
RHS3 = NaN(n_h,1);

for h = 1:n_h
    LHS3(h) = (h-1)^2 * (sigma^2 * rho^(2*(h-1)) - 2 * sigma^2 * rho^(2*h) + rho^(2*(h+1)) * sigma^2 + rho^(2*(h-1)) - 2 * rho^(2*h) + rho^(2*(h+1)));
    RHS3(h) = (1-rho^(2*(h-1))) * (sigma^4 + sigma^2 * rho^2 + sigma^2 + rho^2);
end

min(LHS3(2:end) < RHS3(2:end))

LHS4 = NaN(n_h,1);
RHS4 = NaN(n_h,1);

for h = 1:n_h
    LHS4(h) = (h-1)^2 * rho^(2*(h-1)) * (1-rho^2)^2 * (1 + sigma^2);
    RHS4(h) = (1-rho^(2*(h-1))) * (sigma^4 + sigma^2 * rho^2 + sigma^2 + rho^2);
end

min(LHS4(2:end) < RHS4(2:end))