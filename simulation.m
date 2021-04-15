% Generate random test data
y = test_data();

% Hyperparameters of the MCMC algorithm
T = 100000;                 % the total number of time steps

prob_birth = 0.025;         % probability to do a birth step
sigma_sigma_eps = 0.05;     % the std of the std of the residuals

prob_real = 0.5;            % probability that a root is real
sigmaz = 1;                 % the std of z (in the generation of roots)
pmax = 4;                   % the maximum model order

% Run the MCMC algorithm
[roots_log, p_log, sigma_eps_log] = main(y, T, prob_birth, sigma_sigma_eps, prob_real, sigmaz, pmax, "test");