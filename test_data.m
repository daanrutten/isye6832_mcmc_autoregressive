function [y, roots, sigma_eps, p] = test_data(n, sigma_sigma_eps, sigmaz, pmax)

roots = zeros(pmax, 1);
sigma_eps = abs(normrnd(0, sigma_sigma_eps));
p = randsample(pmax, 1);

for k = 1:floor(0.5*p)
    if rand() < 0.5
        % Sample real roots
        z1z2 = normrnd(0, sigmaz, 2, 1);
        
        % Compute the new roots
        roots(2*k-1:2*k) = 2 * exp(z1z2) ./ (1 + exp(z1z2)) - 1;
    else
        % Sample complex roots
        z1 = abs(normrnd(0, sigmaz));
        u1 = 2 * pi * rand();
        
        % Compute the new roots
        roots(2*k-1:2*k) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * [u1; -u1]);
    end
end

if mod(p, 2) == 1
    % Sample real root
    z = normrnd(0, sigmaz);
    
    % Compute the new root
    roots(p) = 2 * exp(z) ./ (1 + exp(z)) - 1;
end

% Generate residuals
eps = normrnd(0, sigma_eps, n, 1);
y = compute_y(eps, roots, p);