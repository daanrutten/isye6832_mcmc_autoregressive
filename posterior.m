function prob = posterior(roots, roots_z, sigma_eps, p)

eps = compute_eps(y, roots, p);
prior_roots = zeros(0.5 * p, 1);

for k = 1:0.5 * p
    if imag(roots(2*k-1)) == 0
        prior_roots(k) = prob_real * prod(normpdf(roots_z(2*k-1:2*k), 0, sigmaz));
    else
        prior_roots(k) = 2 * (1 - prob_real) * prod(normpdf(roots_z(2*k-1), 0, sigmaz));
    end
end

prob = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * prod(prior_roots) * sigma_eps^(-2);