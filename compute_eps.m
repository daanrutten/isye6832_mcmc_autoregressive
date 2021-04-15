function eps = compute_eps(y, roots, p, pmax)

eps = y;

for k = 1:p
    bin = nchoosek(roots, k);
    eps(pmax+1:end) = eps(pmax+1:end) + sum(prod(-bin, 2), 1) * y(pmax+1-k:end-k);
end