% Generate random test data
y = test_data();

% Hyperparameters of MCMC algorithm
prob_update = 0.95;     % probability to do an update step
prob_birth = 0.025;     % probability to do a birth step
prob_death = 0.025;     % probability to do a death step

T = 100000;             % the total number of time steps
prob_real = 0.5;        % probability that a root is real
sigmaz = 1;             % the std of z (in the generation of roots)
pmax = 4;

% Initial values for parameters
roots = zeros(pmax, 1);     % the current roots (at time t)
prior_roots = ones(pmax, 1);    % prior of the current roots (at time t)
drawn_roots = ones(pmax, 1);    % q of the current roots (at time t)
p = 2;                      % the current p (at time t)

eps = compute_eps(y, roots, p);     % the current residuals (at time t)
sigma_eps = 0.05;                   % the current std of the residuals (at time t)

% Keep values of parameters over time
roots_log = zeros(T, pmax);
p_log = zeros(T, 1);
sigma_eps_log = zeros(T, 1);

for t = 1:T
    probm = compute_probm(prob_update, prob_birth, prob_death, p, pmax);
    u = rand();
    
    if u <= probm(1)
        % Do an update step
        disp("Doing an update step");
        k = randsample(0.5 * p, 1);     % sample a pair of roots to change
        
        % Hyperparameters (TODO: change dynamically)
        mur = 0;
        sigmar = sigmaz;
        
        if rand() < prob_real
            % Sample random z1 and z2 according to normal distribution
            z1 = normrnd(mur, sigmar);
            z2 = normrnd(mur, sigmar);

            % Compute the new roots
            roots_star = roots;
            roots_star(2*k-1:2*k) = 2 * exp([z1; z2]) ./ (1 + exp([z1; z2])) - 1;
            prior_roots_star = prior_roots;
            prior_roots_star(2*k-1:2*k) = normpdf([z1; z2], 0, sigmaz);
            prior_roots_star(2*k-1) = prob_real * prior_roots_star(2*k-1);
            drawn_roots_star = drawn_roots;
            drawn_roots_star(2*k-1:2*k) = normpdf([z1; z2], mur, sigmar);
            drawn_roots_star(2*k-1) = prob_real * drawn_roots_star(2*k-1);
        else
            % Sample random z1 according to normal distribution and u1 uniform
            z1 = abs(normrnd(mur, sigmar));
            u1 = rand();

            % Compute the new roots
            roots_star = roots;
            roots_star(2*k-1:2*k) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * 2 * pi * [u1; -u1]);
            prior_roots_star = prior_roots;
            prior_roots_star(2*k-1) = 2 * (1 - prob_real) * normpdf(z1, 0, sigmaz);
            prior_roots_star(2*k) = 1;
            drawn_roots_star = drawn_roots;
            drawn_roots_star(2*k-1) = 2 * (1 - prob_real) * normpdf(z1, mur, sigmar);
            drawn_roots_star(2*k) = 1;
        end
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p);
        
        % Compute the probabilities used in computing the acceptance ratio
        % (TODO: optimize away constants after verification of correctness)
        pi_mstar_to_m = probm(1) / (0.5 * p) * prod(drawn_roots(2*k-1:2*k));
        pi_m_to_mstar = probm(1) / (0.5 * p) * prod(drawn_roots_star(2*k-1:2*k));        
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps)) * prod(prior_roots_star(2*k-1:2*k));
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * prod(prior_roots(2*k-1:2*k));
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        disp("accratio is " + accratio);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            prior_roots = prior_roots_star;
            drawn_roots = drawn_roots_star;
            eps = eps_star;
            
            disp("Accepted");
        else
            disp("Rejected");
        end
        
        % Sample random sigma according to log normal distribution
        sigma_eps_star = lognrnd(0, 1);
        
        % Compute the probabilities used in computing the acceptance ratio
        % (TODO: optimize away constants after verification of correctness)
        pi_mstar_to_m = probm(1) * lognpdf(sigma_eps, 0, 1);
        pi_m_to_mstar = probm(1) * lognpdf(sigma_eps_star, 0, 1);
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps_star)) * sigma_eps_star^(-2);
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * sigma_eps^(-2);
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        disp("accratio is " + accratio);
        
        if rand() < accratio
            % Move to the new state
            sigma_eps = sigma_eps_star;
            
            disp("Accepted");
        else
            disp("Rejected");
        end
    elseif u <= probm(1) + probm(2)
        % Do a birth step
        disp("Doing a birth step");
        
        % Hyperparameters (TODO: change dynamically) 
        mur = 0;
        sigmar = sigmaz;
        
        if rand() < prob_real
            % Sample random z1 and z2 according to normal distribution
            z1 = normrnd(mur, sigmar);
            z2 = normrnd(mur, sigmar);

            % Compute the new roots
            roots_star = roots;
            roots_star(p+1:p+2) = 2 * exp([z1; z2]) ./ (1 + exp([z1; z2])) - 1;
            prior_roots_star = prior_roots;
            prior_roots_star(p+1:p+2) = normpdf([z1; z2], 0, sigmaz);
            prior_roots_star(p+1) = prob_real * prior_roots_star(p+1);
            drawn_roots_star = drawn_roots;
            drawn_roots_star(p+1:p+2) = normpdf([z1; z2], mur, sigmar);
            drawn_roots_star(p+1) = prob_real * drawn_roots_star(p+1);
        else
            % Sample random z1 according to normal distribution and u1
            % uniform
            z1 = abs(normrnd(mur, sigmar));
            u1 = rand();

            % Compute the new roots
            roots_star = roots;
            roots_star(p+1:p+2) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * 2 * pi * [u1; -u1]);
            prior_roots_star = prior_roots;
            prior_roots_star(p+1) = 2 * (1 - prob_real) * normpdf(z1, 0, sigmaz);
            prior_roots_star(p+2) = 1;
            drawn_roots_star = drawn_roots;
            drawn_roots_star(p+1) = 2 * (1 - prob_real) * normpdf(z1, mur, sigmar);
            drawn_roots_star(p+2) = 1;
        end
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p+2);
        
        % Compute the move probabilities for the new p
        probm_star = compute_probm(prob_update, prob_birth, prob_death, p+2, pmax);
        
        % Compute the probabilities used in computing the acceptance ratio
        % (TODO: optimize away constants after verification of correctness)
        pi_mstar_to_m = probm_star(3) / (0.5 * p + 1);
        pi_m_to_mstar = probm(2) * prod(drawn_roots_star(p+1:p+2));
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps)) * prod(prior_roots_star(p+1:p+2));
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps));
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        disp("accratio is " + accratio);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            prior_roots = prior_roots_star;
            drawn_roots = drawn_roots_star;
            eps = eps_star;
            p = p + 2;
            
            disp("Accepted");
        else
            disp("Rejected");
        end
    else
        % Do a death step
        disp("Doing a death step");
        k = randsample(0.5 * p, 1);     % sample a pair of roots to remove
        
        % Hyperparameters (TODO: change dynamically) 
        mur = 0;
        sigmar = sigmaz;
        
        % Swap out the selected root to the end
        swap = [1:2*(k-1), 2*k+1:pmax, 2*k-1, 2*k]';
        
        roots_star = roots(swap);
        prior_roots_star = prior_roots_star(swap);
        prior_roots_star(pmax-1:pmax) = 1;
        drawn_roots_star = drawn_roots_star(swap);
        drawn_roots_star(pmax-1:pmax) = 1;
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p-2);
        
        % Compute the move probabilities for the new p
        probm_star = compute_probm(prob_update, prob_birth, prob_death, p-2, pmax);
        
        % Compute the probabilities used in computing the acceptance ratio
        % (TODO: optimize away constants after verification of correctness)
        pi_mstar_to_m = probm_star(2) * prod(drawn_roots(2*k-1:2*k));
        pi_m_to_mstar = probm(3) / (0.5 * p);
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps));
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * prod(prior_roots(2*k-1:2*k));
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        disp("accratio is " + accratio);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            prior_roots = prior_roots_star;
            drawn_roots = drawn_roots_star;
            eps = eps_star;
            p = p - 2;
            
            disp("Accepted");
        else
            disp("Rejected");
        end
    end
    
    roots_log(t, :) = roots;
    p_log(t) = p;
    sigma_eps_log(t) = sigma_eps;
end