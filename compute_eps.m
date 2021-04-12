function eps = compute_eps(y, roots, p)

eps = y;

for k = 1:p
    bin = nchoosek(roots, k);
    eps = eps + sum(prod(-bin, 2), 1) * [repelem(0, k)'; y(1:end-k)];
end