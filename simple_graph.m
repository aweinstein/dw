function [L,T] = simple_graph()
% function [L,T] = ring_graph(N)
%
% Return the Laplacian L and the symmetric diffusion operator of a graph with
% 8 vertices that looks like
%
% 1 -- 2 -- 5 -- 6
% |    |    |    |
% |    |    |    |
% 4 -- 3    8 -- 7


W = [0 1 0 1 0 0 0 0;
     1 0 1 0 1 0 0 0;
     0 1 0 1 0 0 0 0;
     1 0 1 0 0 0 0 0;
     0 1 0 0 0 1 0 1;
     0 0 0 0 1 0 1 0;
     0 0 0 0 0 1 0 1;
     0 0 0 0 1 0 1 0];

D = diag(sum(W,2));
 
T = sqrtm(D)^-1 * W * sqrtm(D)^-1;
L = D - W;
