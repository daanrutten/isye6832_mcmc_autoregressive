function [roots_log, p_log, sigma_eps_log] = main(y, T, prob_birth, sigma_sigma_eps, prob_real, sigmaz, pmax)

% Set initial values of parameters
roots = zeros(pmax, 1);                         % the current inverse roots (at time t)
sigma_eps = abs(normrnd(0, sigma_sigma_eps));   % the current std of the residuals (at time t)
p = 2 * randsample(0.5 * pmax, 1);              % the current model order (at time t)

roots_z = zeros(pmax, 1);                       % the z of the current roots (at time t)
prior_roots = normpdf(roots_z, 0, sigmaz);      % prior of the current roots (at time t)
q_roots = prior_roots;                          % q of the current roots (at time t)

% Track values of parameters over time
roots_log = zeros(T, pmax);
sigma_eps_log = zeros(T, 1);
p_log = zeros(T, 1);

% Compute the residuals
eps = compute_eps(y, roots, p);                 % the current residuals (at time t)

for t = 1:T
    probm = compute_probm(prob_birth, p, pmax);
    u = rand();
    
    if u <= probm(1)
        % Do an update step
        k = randsample(0.5 * p, 1);             % sample a pair of roots to change
        
        if rand() < prob_real
            % Sample random z1 and z2 according to normal distribution
            z1z2 = normrnd(0, sigmaz, 2, 1);

            % Compute the new roots
            roots_star = roots;
            roots_star(2*k-1:2*k) = 2 * exp(z1z2) ./ (1 + exp(z1z2)) - 1;
            roots_z_star = roots_z;
            roots_z_star(2*k-1:2*k) = z1z2;
            prior_roots_star = prior_roots;
            prior_roots_star(2*k-1:2*k) = normpdf(z1z2, 0, sigmaz);
            prior_roots_star(2*k-1) = prob_real * prior_roots_star(2*k-1);
            q_roots_star = q_roots;
            q_roots_star(2*k-1:2*k) = normpdf(z1z2, 0, sigmaz);
            q_roots_star(2*k-1) = prob_real * q_roots_star(2*k-1);
        else
            % Sample random z1 according to normal distribution and u1 uniform
            z1 = abs(normrnd(0, sigmaz));
            u1 = 2 * pi * rand();

            % Compute the new roots
            roots_star = roots;
            roots_star(2*k-1:2*k) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * [u1; -u1]);
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
        
        % Compute the acceptance ratio
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
        end
        
        % Sample random sigma according to half-normal distribution
        sigma_eps_star = abs(normrnd(0, sigma_sigma_eps));
        
        % Compute the acceptance ratio
        pi_mstar_to_m = normpdf(sigma_eps, 0, sigma_sigma_eps);
        pi_m_to_mstar = normpdf(sigma_eps_star, 0, sigma_sigma_eps);
        pi_y_mstar = prod(normpdf(eps(pmax+1:end), 0, sigma_eps_star)) * sigma_eps_star^(-2);
        pi_y_m = prod(normpdf(eps(pmax+1:end), 0, sigma_eps)) * sigma_eps^(-2);
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        
        if rand() < accratio
            % Move to the new state
            sigma_eps = sigma_eps_star;
        end
    elseif u <= probm(1) + probm(2)
        % Do a birth step
        if rand() < prob_real
            % Hyperparameters
        	%sigmar = 4 * sigma_eps^2 * sigmaz^2 / (4 * sigma_eps^2 + sigmaz^2 * sum(eps(pmax:end-1).^2));
            %mur = sigmar * sum(eps(pmax:end-1) .* eps(pmax+1:end)) / (2 * sigma_eps^2);
            %sigmar = sqrt(sigmar);
            
            sigmar = sigmaz;
            mur = 0;
            
            % Sample random z1 and z2 according to normal distribution
            z1z2 = normrnd(mur, sigmar, 2, 1);
            
            % Compute the new roots
            roots_star = roots;
            roots_star(p+1:p+2) = 2 * exp(z1z2) ./ (1 + exp(z1z2)) - 1;
            roots_z_star = roots_z;
            roots_z_star(p+1:p+2) = z1z2;
            prior_roots_star = prior_roots;
            prior_roots_star(p+1:p+2) = normpdf(z1z2, 0, sigmaz);
            prior_roots_star(p+1) = prob_real * prior_roots_star(p+1);
            q_roots_star = q_roots;
            q_roots_star(p+1:p+2) = normpdf(z1z2, 0, sigmaz);
            q_roots_star(p+1) = prob_real * q_roots_star(p+1);
            
            b_roots_star = prob_real * prod(normpdf(z1z2, mur, sigmar));
        else
            % Hyperparameters
        	%sigmar = 2 * sigma_eps^2 * sigmaz^2 / (2 * sigma_eps^2 + sigmaz^2 * (2 * sum(eps(pmax:end-1).^2) + sum(eps(pmax-1:end-2) .* eps(pmax+1:end))));
            %mur = sigmar * sum(eps(pmax:end-1) .* eps(pmax+1:end)) / sigma_eps^2;
            %sigmar = sqrt(sigmar);
            
            sigmar = sigmaz;
            mur = 0;
            
            % Sample random z1 according to normal distribution and u1
            % uniform
            z1 = abs(normrnd(mur, sigmar));
            u1 = 2 * pi * rand();
            
            % Compute the new roots
            roots_star = roots;
            roots_star(p+1:p+2) = (2 * exp(z1) ./ (1 + exp(z1)) - 1) * exp(1i * [u1; -u1]);
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
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = 1 / (0.5 * p + 1);
        pi_m_to_mstar = b_roots_star;
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
        end
    else
        % Do a death step
        k = randsample(0.5 * p, 1);         % sample a pair of roots to remove
        
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
            %mur = sigmar * sum(eps_star(pmax:end-1) .* eps_star(pmax+1:end)) / (2 * sigma_eps^2);
            %sigmar = sqrt(sigmar);
            
            sigmar = sigmaz;
            mur = 0;
            
            % Compute the new roots
            b_roots_star = prob_real * prod(normpdf(roots_z(2*k-1:2*k), mur, sigmar));
        else
            % Hyperparameters
        	%sigmar = 2 * sigma_eps^2 * sigmaz^2 / (2 * sigma_eps^2 + sigmaz^2 * (2 * sum(eps_star(pmax:end-1).^2) + sum(eps_star(pmax-1:end-2) .* eps_star(pmax+1:end))));
            %mur = sigmar * sum(eps_star(pmax:end-1) .* eps_star(pmax+1:end)) / sigma_eps^2;
            %sigmar = sqrt(sigmar);
            
            sigmar = sigmaz;
            mur = 0;

            % Compute the new roots
            b_roots_star = 2 * (1 - prob_real) * normpdf(roots_z(2*k-1), mur, sigmar);
        end
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = b_roots_star;
        pi_m_to_mstar = 1 / (0.5 * p);
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
        end
    end
    
    roots_log(t, :) = roots;
    sigma_eps_log(t) = sigma_eps;
    p_log(t) = p;
end