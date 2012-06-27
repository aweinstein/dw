% Example used in 
% M. Maggioni and S. Mahadevan, A Multiscale Framework for Markov Decision
% Processes using Diffusion Wavelets. Amherst: , 2006, pp. 1-44.
% and
% M. Maggioni and S. Mahadevan, "Fast direct policy evaluation using
% multiscale analysis of Markov diffusion processes," in Proceedings of the 
% 23rd international conference on Machine learning, 2006, pp. 601–608.

%close all
clear all
clc

T = [0.8 0.2 0 0;
     0.2 0.75 0.05 0;
     0 0.05 0.75 0.2;
     0 0 0.2 0.8];

j_max = 6;
epsilon = 1e-10;
Tree = DWPTree (T, j_max, epsilon, ...
                struct('Wavelets', true, 'OpThreshold', 1e-5, ...
                       'Symm', true)); %, ...      
                       %'GSOptions', struct('StopDensity', 10, 'Threshold', 1e-3)));
                   

figure(1), clf
js = [1 5 6];
subplot(length(js) + 1, 2, 1)
imagesc(T)
title('T_0')
subplot(length(js) + 1, 2, 2 + length(js))
plot(eye(size(T)))
title('phi_0')
for i = 1:length(js),
    j = js(i);
    T_size = size(Tree{j,1}.T{1},1);
    W_size = size(Tree{j,2}.Basis);
    fprintf('j: %d, T size: %d, W size: (%d, %d)\n', ...
            j, T_size(1), W_size(1), W_size(2))
    subplot(length(js) + 1, 2, i + 1)
    imagesc(Tree{j,1}.T{1})
    title(sprintf('T%d', j))
    subplot(length(js) + 1, 2, i + 2 + length(js))
    plot(Tree{j,1}.ExtBasis)
    title(sprintf('phi_%d', j))
end


%%
% figure
% hold on
% colors = 'bgrcmyk';
% for j = 1:j_max,
%     A = T^(2^(j-1));
%     sigma = svd(A);
%     plot(sort(sigma, 'descend'), [colors(j) '--.'])
% end
