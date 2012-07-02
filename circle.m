% Sec 8.1 (p. 46) of the paper

clear all
close all
clc


use_self_loop = true;
if use_self_loop,
    load T_circle_selfloop_256
    fn = 'circle_tree_sl';
else
    load T_circle_256 % Created by make_T_circle.py
    fn = 'circle_tree';
end

%D = diag(sum(W,1));
%T = D^-0.5 * W * D^-0.5;
load_file = false;

if load_file,
    load(fn)
    fprintf('Loaded %s\n', fn)
else
    GSOptions = struct('StopDensity',1,'Threshold',1e-2);
    opts = struct('Wavelets', true, 'OpThreshold', 1e-2, 'GSOptions', GSOptions);
    tic
    Tree = DWPTree(T, 15, 1e-10, GSOptions); 
    toc
    save(fn, 'Tree')
    fprintf('Tree saved in %s\n', fn)
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

%f_hat = best_k_approx(Tree, wc, 2);

%% Plot eigenvalues of T

figure
evals = eig(T);
figure(1)
plot(sort(evals,'descend'))
grid
xlim([1,size(T,1)])
%ylim([0,1])


%% Plot scaling functions
figure
for j = 1:9, %size(Tree,1),
       
    subplot(3,3,j);
    plot(Tree{j,1}.ExtBasis(:,1))
    hold on
    plot(Tree{j,1}.ExtBasis(:,end/2), 'r')
    title(sprintf('Diffusion scaling function at scale %d',j));
end