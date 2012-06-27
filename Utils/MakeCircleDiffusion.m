function T = MakeCircleDiffusion(N)
% function T = MakeCircleDiffusion(N)
%
% MAKECIRCLEDIFFUSION returns the 3-point diffusion operator on the circle.
%
% In:
%    N = number of points for the circle' discretization
%
% Out:
%    T = NxN sparse matrix representing the diffusion
%

e = ones(N,1);

T = spdiags([0.25*e 0.5*e 0.25*e], -1:1,N,N);

return;


