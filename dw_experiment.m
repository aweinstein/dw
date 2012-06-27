% Experiment using Diffusion Wavelets and Spectral Graph Wavelets
clear all
close all
clc

graph_type = 4;
switch graph_type,
    case 1, % Ring graph
        j_max = 10;
        j_plot = [2 3 5 8];
        N = 128;
        [L, T] = ring_graph(N);
    case 2, % Simple graph
        j_max = 4;
        j_plot = 1:j_max;
        [L, T] = simple_graph();
    case 3,  % T used in Maggioni and Mahadevan paper's example
        j_max = 6;
        j_plot = 1:j_max;
        T = [0.8 0.2 0 0;
             0.2 0.75 0.05 0;
             0 0.05 0.75 0.2;
             0 0 0.2 0.8];
        L = T; % I don't have L. Is this Kocher?
    case 4, % two rings of 64 vertices connected at a single point
        load two_rings
        j_max = 10;
        j_plot = [2 3 5 8];
        T = sqrtm(D)^-1 * A * sqrtm(D)^-1;
    case 5, % five rings of 32 vertices connected at a single point
        load five_rings
        j_max = 10;
        j_plot = [2 3 5 8];
        T = sqrtm(D)^-1 * A * sqrtm(D)^-1;
end

%%%%%%%%%%%% Diffusion Wavelets %%%%%%%%%%%%%%%%%%%
Tree = DWPTree (T, j_max, 1e-4, ...
                struct('Wavelets', true, 'OpThreshold', 1e-2, ...
                'GSOptions', struct('StopDensity', 10, 'Threshold', 1e-3), ...
                'Symm', true));      

% Wavelet coefficient for a delta function            
delta = zeros(size(T,1), 1);
delta(1) = 1;
DWC = DWCoeffs(Tree, delta);

%% Display sizes of V_j and W_j
for j = 1:size(Tree,1),
    T_size = size(Tree{j,1}.T{1},1);
    W_size = size(Tree{j,2}.Basis);
    fprintf('j: %d, T size: %d, W size: (%d, %d)\n', ...
            j, T_size(1), W_size(1), W_size(2))
end

%% Plot T_j
do_plot = true;
if do_plot,
    figure 
    for j = 1:j_max,
        subplot(j_max/2, 2, j)
        imagesc(Tree{j,1}.T{1})
        title(sprintf('T%d', j))
    end
end

%% Plot some of the basis
do_plot = true;
if do_plot,
    figure
    n_rows = ceil(numel(j_plot) / 2);
    for i = 1:numel(j_plot),
        j = j_plot(i);
        a = round(size(Tree{j,1}.T{1},1) / 2);
        subplot(n_rows, 2, i)
        plot(Tree{j,1}.ExtBasis(:,a))
        title(sprintf('phi_%d at location %d', j, a))
        xlim([1 size(T,1)])
    end
end


%%
plot_svds = false;
if plot_svds,
    figure
    hold on
    colors = 'bgrcmyk';
    for j = 1:j_max,
        A = T^(2^(j-1));
        sigma = svd(A);
        col = colors(mod(j-1, length(colors)) + 1);
        plot(sort(sigma, 'descend'), [col '--.'])
    end
end

%%%%%%%%%%%% Spectral Graph Wavelets %%%%%%%%%%%%%%%%%%%
%%
l_max = 1.01 * max(abs(eig(L))); % Use sgwt_rough_lmax(L) for large L
Nscales = 4;
[g, ~, t] = sgwt_filter_design(l_max, Nscales);
m = 50;
arange = [0 l_max];
for k= 1:numel(g),
    c{k} = sgwt_cheby_coeff(g{k}, m, m+1, arange);
end

delta = zeros(size(L,1), 1);
v = round(size(L,1));
delta(v) = 1;
wpall = sgwt_cheby_op(delta, L, c, arange);

%% Plot 
figure
n_rows = ceil(numel(wpall) / 2);
for i = 1:numel(wpall),
    subplot(n_rows, 2, i)
    plot(wpall{i})
    if i == 1,
        title('Scaling function')
    else
        title(sprintf('Wavelt function at scale %d', i-1))
    end
    xlim([1 size(T,1)])
end
