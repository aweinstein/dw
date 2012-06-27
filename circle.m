% Sec 8.1 (p. 46) of the paper

clear all
close all
clc

load Tp
%Tp = sparse(Tp);

T = Tp^2;
evals = eig(T);


load_file = true;

if load_file,
    load('circle_tree')
else
    GSOptions = struct('StopDensity',1,'Threshold',1e-2);
    opts = struct('Wavelets', true, 'OpThreshold', 1e-2, 'GSOptions', GSOptions);
    tic
    Tree = DWPTree(T, 15, 1e-4, GSOptions); 
    toc
    save('circle_tree', 'Tree')
end

% Make function as linear combination of two wavelet functions
j1 = 7;
j2 = 9;
f = 0.5*Tree{j1,2}.ExtBasis(:,1) + Tree{j2,2}.ExtBasis(:,end-3);

%% Take DW transform of f
    wc_tree = DWCoeffs(Tree, f);
    wc_ = DWWavelet(wc_tree); % These are the N coefficients: 7 scaling + 505 wavelets
    wc = DWUnpack(wc_);
    wc = full(wc);

    f_hat = best_k_approx(Tree, wc, 2);

%% Plot

figure
plot(f)
% figure(1)
% plot(sort(evals,'descend'))
% grid
% xlim([1,512])
% ylim([0,1])



% figure
% for j = 1:9,
%     subplot(3,3,j);
%     plot(Tree{j,1}.ExtBasis(:,1))
%     hold on
%     title(sprintf('Diffusion scaling function at scale %d',j));
% end