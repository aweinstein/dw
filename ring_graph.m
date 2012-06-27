function [L,T] = ring_graph(N)
% function [L,T] = ring_graph(N)
%
% Return the Laplacian L and the symmetric diffusion operator of a ring graph 
% with N vertices

c = zeros(1, N);
c(2) = 1;
c(N) = 1;
W = circulant(c)';
D = diag(sum(W,2));
 
T = sqrtm(D)^-1 * W * sqrtm(D)^-1;
L = D - W;
