function y = test_data()

eps = normrnd(0, 0.01, 100, 1);
y = ones(100, 1);

for k = 3:100
    y(k) = 0.25 * y(k-2) + eps(k);  % the roots are -0.5 and +0.5
end