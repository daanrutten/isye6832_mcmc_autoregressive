function out = compute_probm(prob_birth, roots, p, pmax)

prob_br = (1 - 0.5 * mod(p, 2)) * prob_birth * (p <= pmax - 1);
prob_bc = prob_birth * (p <= pmax - 2);
prob_d = prob_birth * (p >= 3 || (p == 2 && imag(roots(1)) == 0));
prob_u = 1 - prob_br - prob_bc - prob_d;

out = [prob_u, prob_br, prob_bc, prob_d];