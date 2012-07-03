clear all
%close all
clc
%%
%load Tree_grid_16
load Tree_grid_sl_16
load Vs
Vs1 = Vs(:,3);
Vsm = reshape(Vs1, 16, 16);

% wc_tree = DWCoeffs(Tree, Vs1);
% Basis  = DWBest(wc_tree);
% wc   = DWUnpack(Basis);
% wc_best = full(wc);

wc_tree = DWCoeffs(Tree, Vs1);
wc_ = DWWavelet(wc_tree); 
wc = DWUnpack(wc_);
wc = full(wc);

wc_tree_ = DWUnpack(wc_tree);
wc_tree_ = full(wc_tree_);

coeffs = sort(abs(wc(:,4)), 'descend');

k = 40;
Vs_hat = best_k_approx(Tree, wc, k);
err_dw = norm(Vs1 - Vs_hat)

% Vs_hat_2 = best_k_approx(Tree, wc_best, k);
% err_dw_2 = norm(Vs1 - Vs_hat_2)



%% k-best approx using the eigenvectors of the laplacian
load T_grid_16
[evecs, evals] = eig(L);
A = evecs(:,1:k);
Vhat_lap = A * pinv(A) * Vs1;
err_lap = norm(Vs1 - Vhat_lap)


%%
figure(1), clf, axis square
subplot(121)
%surface(Vsm)
surf(Vsm)
xlim([1,16]), ylim([1,16])
axis ij
axis square
title('Vs')
subplot(122)
%surface(reshape(Vs_hat,16,16))
surf(reshape(Vs_hat,16,16))
xlim([1,16]), ylim([1,16])
axis ij
axis square
title(sprintf('Vs_{hat} k=%d', k))
%colorbar

%%
figure(2), clf
k_max = 50;
subplot(211)
stem(coeffs)
title('Sorted absolute DW coefficients')
xlim([1,k_max])
subplot(212)
h = evecs' * Vs1;
stem(abs(h))
title('Smoothest absolute Laplacian coefficients')
xlim([1,k_max])

%% Plot some basis

% plot_basis = true;
% if plot_basis,
%     figure(3), clf
%     j = 7;
%     k = 1;
%     i = 1;
%     phi = Tree{j,k}.ExtBasis(:,i);
%     surface(reshape(phi,16,16))
%     axis ij
%     axis square
% end