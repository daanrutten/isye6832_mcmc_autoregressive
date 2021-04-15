% Hyperparameters of the MCMC algorithm
T = 10^8;                   % the total number of time steps

prob_birth = 0.025;         % probability to do a birth step
sigma_sigma_eps = 0.01;     % the std of the std of the residuals

prob_real = 0.5;            % probability that a root is real
sigmaz = 10;                % the std of z (in the generation of roots)
pmax = 2;                   % the maximum model order

name = "test";

% Generate random test data
[y, roots, sigma_eps, p] = test_data(100, sigma_sigma_eps, prob_real, sigmaz, pmax);
save(name + "_real", "y", "T", "prob_birth", "sigma_sigma_eps", "prob_real", "sigmaz", "pmax", "roots", "sigma_eps", "p");

% Run the MCMC algorithm
[roots_log, p_log, sigma_eps_log] = main_static(y, T, prob_birth, sigma_sigma_eps, prob_real, sigmaz, pmax, name);