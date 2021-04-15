function y = compute_y(eps, roots, p)

y = eps;
y(1:p) = 1;

for t = p+1:size(eps, 1)
    for k = 1:p
        bin = nchoosek(roots, k);
        y(t) = y(t) - sum(prod(-bin, 2), 1) * y(t-k);
    end
end