function [roots_log, p_log, sigma_eps_log] = main_dynamic(y, T, prob_birth, sigma_sigma_eps, sigmaz, pmax, name)

Tbatch = min(T, 10^7);

% Set initial values of parameters
roots = zeros(pmax, 1);                         % the current inverse roots (at time t)
roots_z = zeros(pmax, 1);                       % the z of the current roots (at time t)
sigma_eps = abs(normrnd(0, sigma_sigma_eps));   % the current std of the residuals (at time t)
p = randsample(pmax, 1);                        % the current model order (at time t)

% Track values of parameters over time
roots_log = zeros(Tbatch, pmax);
sigma_eps_log = zeros(Tbatch, 1);
p_log = zeros(Tbatch, 1);

% Compute the residuals
eps = compute_eps(y, roots, p, pmax);           % the current residuals (at time t)
prior_eps = -sum(eps(pmax+1:end).^2) / (2 * sigma_eps^2);

for t = 1:T
    probm = compute_probm(prob_birth, roots, p, pmax);
    u = rand();
    
    if u <= probm(1)
        % Do an update step
        k = randsample(p, 1);
        
        if imag(roots(k)) == 0
            % Sample real root
            z = normrnd(0, sigmaz);
            
            % Compute the new root
            roots_star = roots;
            roots_star(k) = 2 * exp(z) ./ (1 + exp(z)) - 1;
            roots_z_star = roots_z;
            roots_z_star(k) = z;
        else
            k2 = ceil(0.5 * k);
            
            % Sample complex root
            z = abs(normrnd(0, sigmaz));
            u1 = 2 * pi * rand();
            
            % Compute the new roots
            roots_star = roots;
            roots_star(2*k2-1:2*k2) = (2 * exp(z) ./ (1 + exp(z)) - 1) * exp(1i * [u1; -u1]);
            roots_z_star = roots_z;
            roots_z_star(2*k2-1:2*k2) = [z; u1];
        end
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p, pmax);
        prior_eps_star = -sum(eps_star(pmax+1:end).^2) / (2 * sigma_eps^2);
        
        % Compute the acceptance ratio
        accratio = exp(prior_eps_star - prior_eps);
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            eps = eps_star;
            prior_eps = prior_eps_star;
        end
        
        if p >= 2
            k = randsample(floor(0.5 * p), 1);
            
            if imag(roots(2*k-1)) ~= 0
                % Sample real root
                z1z2 = normrnd(0, sigmaz, 2, 1);
                
                % Compute the new roots
                roots_star = roots;
                roots_star(2*k-1:2*k) = 2 * exp(z1z2) ./ (1 + exp(z1z2)) - 1;
                roots_z_star = roots_z;
                roots_z_star(2*k-1:2*k) = z1z2;
            else
                % Sample complex root
                z = abs(normrnd(0, sigmaz));
                u1 = 2 * pi * rand();
                
                % Compute the new roots
                roots_star = roots;
                roots_star(2*k-1:2*k) = (2 * exp(z) ./ (1 + exp(z)) - 1) * exp(1i * [u1; -u1]);
                roots_z_star = roots_z;
                roots_z_star(2*k-1:2*k) = [z; u1];
            end
            
            % Compute the residuals for the new roots
            eps_star = compute_eps(y, roots_star, p, pmax);
            prior_eps_star = -sum(eps_star(pmax+1:end).^2) / (2 * sigma_eps^2);
            
            % Compute the acceptance ratio
            accratio = exp(prior_eps_star - prior_eps);
            
            if rand() < accratio
                % Move to the new state
                roots = roots_star;
                roots_z = roots_z_star;
                eps = eps_star;
                prior_eps = prior_eps_star;
            end
        end
        
        % Sample sigma
        sigma_eps_star = abs(normrnd(0, sigma_sigma_eps));
        prior_eps_star = -sum(eps(pmax+1:end).^2) / (2 * sigma_eps_star^2);
        
        % Compute the acceptance ratio
        pi_mstar_to_m = normpdf(sigma_eps, 0, sigma_sigma_eps);
        pi_m_to_mstar = normpdf(sigma_eps_star, 0, sigma_sigma_eps);
        pi_y_mstar = (sigma_eps / sigma_eps_star)^(size(y, 1) - pmax) * exp(prior_eps_star - prior_eps) * sigma_eps_star^(-2);
        pi_y_m = sigma_eps^(-2);
        
        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);
        
        if rand() < accratio
            % Move to the new state
            sigma_eps = sigma_eps_star;
            prior_eps = prior_eps_star;
        end
    elseif u <= probm(1) + probm(2)
        % Do a real birth step
        
        % Hyperparameters
        sigmar = 4 * sigma_eps^2 * sigmaz^2 / (4 * sigma_eps^2 + sigmaz^2 * sum(eps(pmax:end-1).^2));
        mur = sigmar * sum(eps(pmax:end-1) .* eps(pmax+1:end)) / (2 * sigma_eps^2);
        sigmar = sqrt(sigmar);
        
        % Sample real root
        z = normrnd(mur, sigmar);
        
        % Compute the new root
        roots_star = roots;
        roots_star(p+1) = 2 * exp(z) ./ (1 + exp(z)) - 1;
        roots_z_star = roots_z;
        roots_z_star(p+1) = z;
        p_star = p+1;
        
        prior_roots_star = normpdf(z, 0, sigmaz);
        q_roots_star = normpdf(z, mur, sigmar);
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p_star, pmax);
        prior_eps_star = -sum(eps_star(pmax+1:end).^2) / (2 * sigma_eps^2);
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = 1 / p_star;
        pi_m_to_mstar = q_roots_star;
        pi_y_mstar = exp(prior_eps_star - prior_eps) * prior_roots_star;
        
        accratio = pi_y_mstar * pi_mstar_to_m / pi_m_to_mstar;
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            eps = eps_star;
            prior_eps = prior_eps_star;
            p = p_star;
        end
    elseif u <= probm(1) + probm(2) + probm(3)
        % Do a complex birth step
        
        % Hyperparameters
        sigmar = 2 * sigma_eps^2 * sigmaz^2 / (2 * sigma_eps^2 + sigmaz^2 * (2 * sum(eps(pmax:end-1).^2) + sum(eps(pmax-1:end-2) .* eps(pmax+1:end))));
        mur = sigmar * sum(eps(pmax:end-1) .* eps(pmax+1:end)) / sigma_eps^2;
        sigmar = sqrt(sigmar);
        
        % Sample complex root
        z = abs(normrnd(mur, sigmar));
        u1 = 2 * pi * rand();
    
        % Compute the new roots
        p2 = ceil(0.5 * p);
        
        % Swap out the last pair of roots to the end
        swap = [1:2*(p2-1), 2*(p2+1)-1, 2*(p2+1), 2*p2-1, 2*p2, 2*(p2+2)-1:pmax]';
        
        roots_star = roots(swap);
        roots_star(2*p2-1:2*p2) = (2 * exp(z) ./ (1 + exp(z)) - 1) * exp(1i * [u1; -u1]);
        roots_z_star = roots_z(swap);
        roots_z_star(2*p2-1:2*p2) = [z; u1];
        p_star = p+2;
        
        prior_roots_star = normpdf(z, 0, sigmaz);
        q_roots_star = normpdf(z, mur, sigmar);
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p_star, pmax);
        prior_eps_star = -sum(eps_star(pmax+1:end).^2) / (2 * sigma_eps^2);
        
        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = 1 / p_star;
        pi_m_to_mstar = q_roots_star;
        pi_y_mstar = exp(prior_eps_star - prior_eps) * prior_roots_star;
        
        accratio = pi_y_mstar * pi_mstar_to_m / pi_m_to_mstar;
        
        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            eps = eps_star;
            prior_eps = prior_eps_star;
            p = p_star;
        end
    else
        % Do a death step
        k = randsample(p, 1);
        
        if imag(roots(k)) == 0
            if k ~= p
                k2 = ceil(0.5 * k);
                
                % Swap out the selected roots to the end
                swap = [1:2*(k2-1), 2*k2+1:p, 2*k2-1 + mod(k, 2), k, p+1:pmax]';
                roots_star = roots(swap);
                roots_z_star = roots_z(swap);
            else
                roots_star = roots;
                roots_z_star = roots_z;
            end
            
            p_star = p-1;
        else
            k2 = ceil(0.5 * k);
            
            % Swap out the selected roots to the end
            swap = [1:2*(k2-1), 2*k2+1:pmax, 2*k2-1, 2*k2]';
            
            roots_star = roots(swap);
            roots_z_star = roots_z(swap);
            p_star = p-2;
        end
        
        % Compute the residuals for the new roots
        eps_star = compute_eps(y, roots_star, p_star, pmax);
        prior_eps_star = -sum(eps_star(pmax+1:end).^2) / (2 * sigma_eps^2);

        if imag(roots(k)) == 0
            % Hyperparameters
            sigmar = 4 * sigma_eps^2 * sigmaz^2 / (4 * sigma_eps^2 + sigmaz^2 * sum(eps_star(pmax:end-1).^2));
            mur = sigmar * sum(eps_star(pmax:end-1) .* eps_star(pmax+1:end)) / (2 * sigma_eps^2);
            sigmar = sqrt(sigmar);
            
            % Compute the new roots
            prior_roots = normpdf(roots_z(k), 0, sigmaz);
            q_roots = normpdf(roots_z(k), mur, sigmar);
        else
            % Hyperparameters
            sigmar = 2 * sigma_eps^2 * sigmaz^2 / (2 * sigma_eps^2 + sigmaz^2 * (2 * sum(eps_star(pmax:end-1).^2) + sum(eps_star(pmax-1:end-2) .* eps_star(pmax+1:end))));
            mur = sigmar * sum(eps_star(pmax:end-1) .* eps_star(pmax+1:end)) / sigma_eps^2;
            sigmar = sqrt(sigmar);

            % Compute the new roots
            prior_roots = normpdf(roots_z(2*k2-1), 0, sigmaz);
            q_roots = normpdf(roots_z(2*k2-1), mur, sigmar);
        end

        % Compute the probabilities used in computing the acceptance ratio
        pi_mstar_to_m = q_roots;
        pi_m_to_mstar = 1 / p;
        pi_y_mstar = exp(prior_eps_star - prior_eps);
        pi_y_m = prior_roots;

        accratio = pi_y_mstar * pi_mstar_to_m / (pi_y_m * pi_m_to_mstar);

        if rand() < accratio
            % Move to the new state
            roots = roots_star;
            roots_z = roots_z_star;
            eps = eps_star;
            prior_eps = prior_eps_star;
            p = p_star;
        end
    end
    
    step = mod(t - 1, Tbatch) + 1;
    roots_log(step, :) = roots;
    sigma_eps_log(step) = sigma_eps;
    p_log(step) = p;
    
    if mod(t, Tbatch) == 0
        save(name + "_" + (t / Tbatch), "y", "T", "prob_birth", "sigma_sigma_eps", "sigmaz", "pmax", "roots_log", "sigma_eps_log", "p_log");
        disp("Finished running " + t + " timesteps");
    end
end