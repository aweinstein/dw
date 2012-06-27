function errors = Vs_error()

load Vs

N = size(Vs, 1);
errors = zeros(N,1);

for k = 1:2,
    Vs_hat = grid_best_approx(Vs, k);
    errors(k) = norm(Vs - Vs_hat, 'fro');
end
