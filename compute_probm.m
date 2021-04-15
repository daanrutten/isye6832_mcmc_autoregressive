function out = compute_probm(prob_birth, p, pmax)

prob_b = prob_birth * (p < pmax - 1);
prob_d = prob_birth * (p > 2);
prob_u = 1 - prob_b - prob_d;

out = [prob_u, prob_b, prob_d];