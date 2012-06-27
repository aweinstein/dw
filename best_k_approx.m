function f_hat = best_k_approx(Tree, wc, l)

coeffs = wc(:,4);
[~, i_sort_coeffs] = sort(abs(coeffs), 'descend');

N = size(wc, 1);
f_hat = zeros(N, 1);


for i = 1:l,
    wc(i_sort_coeffs(i), :);
    j = wc(i_sort_coeffs(i), 1);
    k = wc(i_sort_coeffs(i), 2);
    col = wc(i_sort_coeffs(i), 3);
    coeff = wc(i_sort_coeffs(i), 4);
    f_hat = f_hat + coeff * Tree{j, k}.ExtBasis(:,col);
end