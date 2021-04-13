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
roots = zeros(pmax, 1);     % the current inverse roots (at time t)
roots_z = zeros(pmax, 1);   % the z of the current roots (at time t)
prior_roots = ones(pmax, 1);    % prior of the current roots (at time t)
q_roots = ones(pmax, 1);        % q of the current roots (at time t)
p = 2;                      % the current p (at time t)

eps = compute_eps(y, roots, p);     % the current residuals (at time t)
sigma_eps = 0.05;                   % the current std of the residuals (at time t)

% Keep values of parameters over time
roots_log = zeros(T, pmax);
p_log = zeros(T, 1);
sigma_eps_log = zeros(T, 1);

acc_update = 0;
acc_update_sigma = 0;
acc_birth = 0;
acc_death = 0;

for t = 1:T
    probm = compute_probm(prob_update, prob_birth, prob_death, p, pmax);
    u = rand();
    
    if u <= probm(1)
        % Do an update step
        k = randsample(0.5 * p, 1);     % sample a pair of roots to change
        
        if rand() < prob_real
            % Sample random z1 and z2 according to normal distribution
            z1 = normrnd(0, sigmaz);
            z2 = normrnd(0, sigmaz);

            % Compute the new roots
            roots_star = roots;
            roots_star(2*k-1:2*k) = 2 * exp([z1; z2]) ./ (1 + exp([z1; z2])) - 1;
            roots_z_star = roots_z;
            roots_z_star(2*k-1:2*k) = [z1; z2];
            prior_roots_star = prior_roots;
            prior_roots_star(2*k-1:2*k) = normpdf([z1; z2], 0, sigmaz);
            prior_roots_star(2*k-1) = prob_real * prior_roots_star(2*k-1);
            q_roots_star = q_roots;
            q_roots_star(2*k-1:2*k) = normpdf([z1; z2], 0, sigmaz);
            q_roots_star(2*k-1) = prob_real * q_roots_star(2*k-1);
        else
            % Sample random z1 according to normal distribution and u1 uniform
            z1 = abs(normrnd(0, sigmaz));
            u1 = rand();

            % Compute the new roots
            roots_star = roots;
            roots_star(2*k-1:2*k) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * 2 * pi * [u1; -u1]);
            roots_z_star = roots_z;
            roots_z_star(2*k-1:2*k) = [z1; u1];
            prior_roots_star = prior_roots;
            prior_roots_star(2*k-1) = 2 * (1 - prob_real) * normpdf(z1, 0, sigmaz);
            prior_roots_star(2*k) = 1;
            q_roots_star = q_roots;
            q_roots_star(2*k-1) = 2 * (1 - prob_real) * normpdf(z1, 0, sigmaz);
            q_roots_star(2*k) = 1;
        end
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p);
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = prod(q_roots(2*k-1:2*k));
        pi_m_to_mstar = prod(q_roots_star(2*k-1:2*k));        
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps)) * prod(prior_roots_star(2*k-1:2*k));
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * prod(prior_roots(2*k-1:2*k));
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            prior_roots = prior_roots_star;
            q_roots = q_roots_star;
            eps = eps_star;
            
            acc_update = acc_update + 1;
        end
        
        % Sample random sigma according to log normal distribution
        sigma_eps_star = abs(normrnd(0, 0.05));
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = normpdf(sigma_eps, 0, 0.05);
        pi_m_to_mstar = normpdf(sigma_eps_star, 0, 0.05);
        pi_y_mstar = prod(normpdf(eps(pmax+1:end), 0, sigma_eps_star)) * sigma_eps_star^(-2);
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * sigma_eps^(-2);
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        
        if rand() < accratio
            % Move to the new state
            sigma_eps = sigma_eps_star;
            
            acc_update_sigma = acc_update_sigma + 1;
        end
    elseif u <= probm(1) + probm(2)
        % Do a birth step
        if rand() < prob_real
            % Hyperparameters
        	%sigmar = 4 * sigma_eps^2 * sigmaz^2 / (4 * sigma_eps^2 + sigmaz^2 * sum(eps(pmax:end-1).^2));
            %mur = sigmar^2 * sum(eps(pmax:end-1) .* eps(pmax+1:end)) / (2 * sigma_eps^2);
            
            sigmar = sigmaz;
            mur = 0;
            
            % Sample random z1 and z2 according to normal distribution
            z1 = normrnd(mur, sigmar);
            z2 = normrnd(mur, sigmar);

            % Compute the new roots
            roots_star = roots;
            roots_star(p+1:p+2) = 2 * exp([z1; z2]) ./ (1 + exp([z1; z2])) - 1;
            roots_z_star = roots_z;
            roots_z_star(p+1:p+2) = [z1; z2];
            prior_roots_star = prior_roots;
            prior_roots_star(p+1:p+2) = normpdf([z1; z2], 0, sigmaz);
            prior_roots_star(p+1) = prob_real * prior_roots_star(p+1);
            q_roots_star = q_roots;
            q_roots_star(p+1:p+2) = normpdf([z1; z2], 0, sigmaz);
            q_roots_star(p+1) = prob_real * q_roots_star(p+2);
            b_roots_star = prob_real * prod(normpdf([z1; z2], mur, sigmar));
        else
            % Hyperparameters
        	%sigmar = 2 * sigma_eps^2 * sigmaz^2 / (2 * sigma_eps^2 + sigmaz^2 * (2 * sum(eps(pmax:end-1).^2) + sum(eps(pmax-1:end-2) .* eps(pmax+1:end))));
            %mur = sigmar^2 * sum(eps(pmax:end-1) .* eps(pmax+1:end)) / sigma_eps^2;
            
            sigmar = sigmaz;
            mur = 0;
            
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
            q_roots_star = q_roots;
            q_roots_star(p+1) = 2 * (1 - prob_real) * normpdf(z1, 0, sigmaz);
            q_roots_star(p+2) = 1;
            b_roots_star = 2 * (1 - prob_real) * normpdf(z1, mur, sigmar);
        end
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p+2);
        
        % Compute the move probabilities for the new p
        probm_star = compute_probm(prob_update, prob_birth, prob_death, p+2, pmax);
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = probm_star(3) / (0.5 * p + 1);
        pi_m_to_mstar = probm(2) * b_roots_star;
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps)) * prod(prior_roots_star(p+1:p+2));
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps));
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            prior_roots = prior_roots_star;
            q_roots = q_roots_star;
            eps = eps_star;
            p = p + 2;
            
            acc_birth = acc_birth + 1;
        end
    else
        % Do a death step
        k = randsample(0.5 * p, 1);     % sample a pair of roots to remove
        
        % Swap out the selected root to the end
        swap = [1:2*(k-1), 2*k+1:pmax, 2*k-1, 2*k]';
        
        roots_star = roots(swap);
        roots_z_star = roots_z(swap);
        prior_roots_star = prior_roots_star(swap);
        q_roots_star = q_roots(swap);
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p-2);
        
        % Compute the q for the new p
        if imag(roots(2*k-1)) == 0
            % Hyperparameters
        	%sigmar = 4 * sigma_eps^2 * sigmaz^2 / (4 * sigma_eps^2 + sigmaz^2 * sum(eps_star(pmax:end-1).^2));
            %mur = sigmar^2 * sum(eps_star(pmax:end-1) .* eps_star(pmax+1:end)) / (2 * sigma_eps^2);
            
            sigmar = sigmaz;
            mur = 0;
            
            % Compute the new roots
            b_roots_star = prob_real * prod(normpdf(roots_z(2*k-1:2*k), mur, sigmar));
        else
            % Hyperparameters
        	%sigmar = 2 * sigma_eps^2 * sigmaz^2 / (2 * sigma_eps^2 + sigmaz^2 * (2 * sum(eps_star(pmax:end-1).^2) + sum(eps_star(pmax-1:end-2) .* eps_star(pmax+1:end))));
            %mur = sigmar^2 * sum(eps_star(pmax:end-1) .* eps_star(pmax+1:end)) / sigma_eps^2;
            
            sigmar = sigmaz;
            mur = 0;

            % Compute the new roots
            b_roots_star = 2 * (1 - prob_real) * normpdf(roots_z(2*k-1), mur, sigmar);
        end
        
        % Compute the move probabilities for the new p
        probm_star = compute_probm(prob_update, prob_birth, prob_death, p-2, pmax);
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = probm_star(2) * b_roots_star;
        pi_m_to_mstar = probm(3) / (0.5 * p);
        pi_y_mstar = prod(normpdf(eps_star(pmax+1:end), 0, sigma_eps));
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * prod(prior_roots(2*k-1:2*k));
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            prior_roots = prior_roots_star;
            q_roots = q_roots_star;
            eps = eps_star;
            p = p - 2;
            
            acc_death = acc_death + 1;
        end
    end
    
    roots_log(t, :) = roots;
    p_log(t) = p;
    sigma_eps_log(t) = sigma_eps;
end