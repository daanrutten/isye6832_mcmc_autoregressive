roots_log_sub = roots_log(p_log == p);
roots_log_sub_real = roots_log_sub(imag(roots_log_sub) == 0);
roots_log_sub_complex = roots_log_sub(imag(roots_log_sub) ~= 0);

figure;
histogram(real(roots_log_sub_real));
title("real roots");

figure;
histogram(real(roots_log_sub_complex));
title("real part of complex roots");

figure;
histogram(imag(roots_log_sub_complex));
title("imaginary part of complex roots");

figure;
histogram(p_log);
title("p");

figure;
histogram(sigma_eps_log);
title("sigma eps");