function F_hat = grid_best_approx(F, k)

% Assume is an M-by-N matrix, where each column of F is a function to be
% approximated

F_hat = zeros(size(F));
N = sqrt(size(F,1));
if ~ (N==4 || N==8 || N==16),
    disp 'Error, cannot handle the dimension'
    return
end

fn = sprintf('Tree_grid_%d', N);
load(fn)

for i = 1:size(F,2),
    wc_tree = DWCoeffs(Tree, F(:,i));
    wc_ = DWWavelet(wc_tree);
    wc = DWUnpack(wc_);
    wc = full(wc);
    F_hat(:,i) = best_k_approx(Tree, wc, k);
end
