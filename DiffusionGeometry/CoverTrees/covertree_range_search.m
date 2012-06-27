function [idxs,distssq] = covertree_range_search( X, CoverTree, Y, R )

%
% Find the R-neighbors of the points in Y in X
%
%
% IN:
%	X			: D by N matrix of N points in R^D
%	CoverTree	: covertree constructed by the covertree function
%	Y			: D by M matrix of M query points
%	R			: radius for the search
%
% OUT:
%	idxs		: M cell array with idxs{i} listing the R-neighbors of the i-th query point
%	dists		: M cell array with idxs{i} listing the distances of the R-neighbors of the i-th query point
%
% EXAMPLE:
% X=[0:127]/128;A.vectors=X+0.1*randn(size(X));A.theta=.5;A.maxlevel=int32(8);CoverTree=covertree(A)
% Y=X(:,CoverTree.levels(1:CoverTree.level_offsets(5))+1);
% [nnIdxs,nnDistssq] = covertree_range_search(X,CoverTree,Y,0.1)

idxs = [];
distssq = [];

C.vectors   = Y;
C.within    = R;
C.level     = int32(CoverTree.params(7))+1;

if nargout>1,
    [idxs] = findwithin(X,CoverTree,C);
    for k = 1:length(idxs),
        distssq_tmp = sum(bsxfun(@minus,X(:,idxs{k}+1),Y(:,k)).^2,1);
        [distssq{k},sorting_idxs] = sort(distssq_tmp,'ascend');
        idxs{k} = idxs{k}(sorting_idxs)+1;
    end;
else
    idxs = findwithin(X,CoverTree,C);
    for k = 1:length(idxs),
        idxs{k} = idxs{k}+1;
    end;    
end;

return;