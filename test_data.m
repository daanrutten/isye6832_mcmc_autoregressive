function [y, roots, sigma_eps, p] = test_data(n, sigma_sigma_eps, prob_real, sigmaz, pmax)

roots = zeros(pmax, 1);
sigma_eps = abs(normrnd(0, sigma_sigma_eps));
p = 2 * randsample(0.5 * pmax, 1);

for k = 1:0.5*p
    if rand() < prob_real
        % Sample random z1 and z2 according to normal distribution
        z1z2 = normrnd(0, sigmaz, 2, 1);
        
        % Compute the new roots
        roots(2*k-1:2*k) = 2 * exp(z1z2) ./ (1 + exp(z1z2)) - 1;
    else
        % Sample random z1 according to normal distribution and u1 uniform
        z1 = abs(normrnd(0, sigmaz));
        u1 = 2 * pi * rand();
        
        % Compute the new roots
        roots(2*k-1:2*k) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * [u1; -u1]);
    end
end

% Generate residuals
eps = normrnd(0, sigma_eps, n, 1);
y = compute_y(eps, roots, p);