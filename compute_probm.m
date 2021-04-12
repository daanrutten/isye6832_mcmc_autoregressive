function y = compute_probm(prob_update, prob_birth, prob_death, p, pmax)

prob_u = prob_update;
prob_b = prob_birth * (p < pmax);
prob_d = prob_death * (p > 2);

y = [prob_u; prob_b; prob_d] / (prob_u + prob_b + prob_d);