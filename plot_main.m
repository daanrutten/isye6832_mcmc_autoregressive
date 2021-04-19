names = ["sim_dym_pb25_pr5_pmax" + (2:2:10), ...
    "sim_dym_pb1_pr5_pmax6", "sim_dym_pb5_pr5_pmax6", ...
    "sim_dym_pb25_pr35_pmax6", "sim_dym_pb25_pr65_pmax6", ...
    "sim_stat_pb25_pr5_pmax" + (2:2:10), ...
    "sim_stat_pb1_pr5_pmax6", "sim_stat_pb5_pr5_pmax6", ...
    "sim_stat_pb25_pr35_pmax6", "sim_stat_pb25_pr65_pmax6"];

for name = names
    % Load variables
    load(name + "_real", "T", "sigma_sigma_eps", "pmax", "roots", "sigma_eps", "p");
    disp("Now processing " + name);
    
    Tbatch = min(T, 10^7);
    
    root_edges = (-1:0.01:1)';
    sigma_edges = (0:0.005*sigma_sigma_eps:2*sigma_sigma_eps)';
    p_edges = (1:pmax+1)';
    
    roots_log_real_total = zeros(1, size(root_edges, 1) - 1);
    roots_log_complexre_total = zeros(1, size(root_edges, 1) - 1);
    roots_log_complexim_total = zeros(1, size(root_edges, 1) - 1);
    sigma_eps_log_total = zeros(1, size(sigma_edges, 1) - 1);
    p_log_total = zeros(1, size(p_edges, 1) - 1);

    for t = 1:20
        if isfile(name + "_" + t + ".mat")
            load(name + "_" + t, "roots_log", "sigma_eps_log", "p_log");
            
            roots_log = roots_log(p_log == p);
            roots_log_real = roots_log(imag(roots_log) == 0);
            roots_log_complex = roots_log(imag(roots_log) ~= 0);
            
            roots_log_real_total = roots_log_real_total + histcounts(real(roots_log_real(:)), root_edges);
            roots_log_complexre_total = roots_log_complexre_total + histcounts(real(roots_log_complex(:)), root_edges);
            roots_log_complexim_total = roots_log_complexim_total + histcounts(imag(roots_log_complex(:)), root_edges);
            sigma_eps_log_total = sigma_eps_log_total + histcounts(sigma_eps_log, sigma_edges);
            p_log_total = p_log_total + histcounts(p_log, p_edges);
        else
            disp("Loaded " + (t-1) + " batches");
            break
        end
    end

    % Plot variables
    set(0, 'defaultfigurecolor', [1 1 1]);
    set(0, 'defaultaxesfontsize', 16);
    set(0, 'defaultaxesticklabelinterpreter', 'latex');
    set(0, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultlegendfontsize', 16);
    set(0, 'defaultlegendinterpreter', 'latex');
    C = linspecer(2);

    figure;
    hold on;
    grid on;
    box on;
    set(gca, "XMinorTick", "on", "YMinorTick", "on");

    bar(root_edges(1:end-1), roots_log_real_total);

    for k = 1:p
        if imag(roots(k)) == 0
            h = xline(real(roots(k)));
            set(h, "Color", C(2, :));
            set(h, "LineWidth", 1.5);
            set(h, "LineStyle", ":");
        end
    end

    xlabel("$\lambda : \lambda \in \bf{R}$");
    ylabel("Frequency");
    xlim([-1 1]);

    savefig(name + "_reallambda");
    saveas(gca, name + "_reallambda.png");

    figure;
    hold on;
    grid on;
    box on;
    set(gca, "XMinorTick", "on", "YMinorTick", "on");

    bar(root_edges(1:end-1), roots_log_complexre_total);

    for k = 1:0.5*p
        if imag(roots(2*k-1)) ~= 0
            h = xline(real(roots(2*k-1)));
            set(h, "Color", C(2, :));
            set(h, "LineWidth", 1.5);
            set(h, "LineStyle", ":");
        end
    end

    xlabel("Re$(\lambda) : \lambda \in \bf{C} \setminus \bf{R}$");
    ylabel("Frequency");
    xlim([-1 1]);

    savefig(name + "_complexlambdare");
    saveas(gca, name + "_complexlambdare.png");

    figure;
    hold on;
    grid on;
    box on;
    set(gca, "XMinorTick", "on", "YMinorTick", "on");

    bar(root_edges(1:end-1), roots_log_complexim_total);

    for k = 1:p
        if imag(roots(k)) ~= 0
            h = xline(imag(roots(k)));
            set(h, "Color", C(2, :));
            set(h, "LineWidth", 1.5);
            set(h, "LineStyle", ":");
        end
    end

    xlabel("Im$(\lambda) : \lambda \in \bf{C} \setminus \bf{R}$");
    ylabel("Frequency");
    xlim([-1 1]);

    savefig(name + "_complexlambdaim");
    saveas(gca, name + "_complexlambdaim.png");

    figure;
    hold on;
    grid on;
    box on;
    set(gca, "XMinorTick", "on", "YMinorTick", "on");

    bar(p_edges(1:end-1), p_log_total);

    h = xline(p);
    set(h, "Color", C(2, :));
    set(h, "LineWidth", 1.5);
    set(h, "LineStyle", ":");

    xlabel("$k$");
    ylabel("Frequency");
    xlim([1 pmax+1]);

    savefig(name + "_k");
    saveas(gca, name + "_k.png");

    figure;
    hold on;
    grid on;
    box on;
    set(gca, "XMinorTick", "on", "YMinorTick", "on");

    bar(sigma_edges(1:end-1), sigma_eps_log_total);

    h = xline(sigma_eps);
    set(h, "Color", C(2, :));
    set(h, "LineWidth", 1.5);
    set(h, "LineStyle", ":");

    xlabel("$\sigma_\varepsilon$");
    ylabel("Frequency");
    xlim([0 2 * sigma_sigma_eps]);

    savefig(name + "_sigma_eps");
    saveas(gca, name + "_sigma_eps.png");
end