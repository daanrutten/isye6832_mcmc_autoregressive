roots_log_2 = roots_log(p_log == 2);
roots_log_2_real = roots_log_2(imag(roots_log_2) == 0);
roots_log_2_complex = roots_log_2(imag(roots_log_2) ~= 0);

figure;
histogram(real(roots_log_2_real));
title("real roots");

figure;
histogram(real(roots_log_2_complex));
title("complex roots");

figure;
histogram(p_log);
title("p");

figure;
histogram(sigma_eps_log);
title("sigma eps");