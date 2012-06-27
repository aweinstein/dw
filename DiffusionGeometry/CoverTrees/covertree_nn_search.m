function [idxs,dists] = covertree_nn_search( X, CoverTree, Y, kNN )

%
% Find the kNN nearest neighbors of the points in Y in X
%
% IN:
%	X			: D by N matrix of N points in R^D
%	CoverTree	: covertree constructed by the covertree function
%	Y			: D by M matrix of M query points
%	kNN			: how many nearest neighbors to compute
%
% OUT:
%	idxs		: M by kNN matrix of indices of nearest neighbors
%	dists		: M by kNN matrix of corresponding distances.
%
%
% EXAMPLE:
% X=[0:127]/128;A.vectors=X+0.1*randn(size(X));A.theta=.5;A.maxdescend=int32(8);CoverTree=covertree(A)
% Y=X(:,CoverTree.levels(1:CoverTree.level_offsets(5))+1);
% [nnIdxs,nnDists] = covertree_nn_search(A.X,CoverTree,Y,2)
%

C.vectors = Y;
C.k       = int32(kNN);

[idxs_cell,dists_cell] = findnearest(X,CoverTree,C);

idxs    = zeros(size(Y,2),kNN);
dists   = zeros(size(Y,2),kNN);
for k = 1:size(idxs,1),
    idxs(k,:)   = idxs_cell{k}(1:kNN)+1;
    dists(k,:)  = dists_cell{k}(1:kNN);
end;

return;