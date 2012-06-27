function [count,idxs,dists,NNInfo] = nrsearch(cX, cXq, NN, radius, metric,opts)

%
% function [count, idxs, dists] = nrsearch(cX, cXq, NN, radius, metric,opts)
%
% NRSEARCH is a wrapper for performing approximate nearest neighbor
% searches using any of the following:
% 1) ANN package of Mount, et al. http://www.cs.umd.edu/~mount/ANN/
% 2) nn_search functions in the OpenTSTool box (http://www.physik3.gwdg.de/tstool/)
% 3) covertrees (to be fully integrated in V2.0 of the Diffusion Geometry toolbox
%
% IN:
%    cX         : DxN matrix of N points in R^D comprising the data set. Could transposed by setting the options XIsTransposed. MUST be of type double.
%                   It may also be an L row vector of indices into cXq if cX is a subset of cXq.
%    cXq        : LxD matrix of L query points in R^D. Must be of type double.
%                   If of type 'uint32', it is a L row vector of indices into cX if cXq is a subset of cX.
%                   This is not an option if D=1, unless DistInfo is provided.
%    NN         : number of nearest neighbors to find
%    radius     : range search with this radius (and at most NN nearest neighbors). Disregarded if 0.
%    metric     : (optional) distance function
%    [opts]     : structure containing the following fields:
%                   [DistRestrict]     : N array of cells, the k-th cell contains the indices of the points among which
%                                        the nearest neighbors of point k are searched for.
%                                        If it is only one cell, then it is assumed to contain the *relative* (w.r.t. current point)
%                                        indices of the points among which the nearest neighbors of point k are searched for.
%                                        Indices out of range are automatically clipped.
%                   [Tolerance]        : tolerance in range search. Default: 0 (exact).
%                   [ReturnDistSquared]: return the distance squared (for Euclidean) distances. Default: 0.
%                   [ReturnAsArrays]   : returns <idxs> and <dists> as arrays rather than as cell arrays. This is only
%                                        possible for NN searches, not for radius searches. Default: 0.
%                   [DistInfo]         : if not empty, it will be used for distance computation. It is assumed that it contains at least
%                                        the neighbors up to the requested distance Radius (and PCARadius if provided), or kNN (if provided).
%                                        It assumes the neighbors are listed in increasing distance. It's a structure containing the following fields:
%                                        In this case cX and cXq are indices, not data points.
%                                        It may be in two formats:
%                                           idxs        : N cell, the k-th cell contains the indices of the nearest neighbors of the k-th point
%                                           dists       : N cell, the k-th cell contains the distances of the nearest neighbors of the k-th point (corresponding to idxs).
%                                        or:
%                                           count       : N vector, k-th entry being number of neighbors of k-th point stored
%                                           idxs        : N by N matrix, the k-th row contaning the indices of neighbors, in order of increasing distance
%                                           dists       : N by N nmatrix, the k-th row containing the distances, in increasing order, from the k-th point
%
%                                     Should be sorted in incrasing magnitude.
%                   [NNInfo]        : extra structure for nearest neighbor extra information.
%                   Atria           : If the TSTOOL nn_prepare is used, this is the atria structure of cX. It will be used for faster computation.
%                   XIsTransposed   : cX is passed transposed, i.e. it is D by N representing N points in D dimensions. Default: 0.
%                   [FastNNSearcher]: 'nn_search','ANNsearch','bruteforce','best'
%                                       (i.e. tries nn_search first, then ANNsearch, then bruteforce).
%                                       Default: 'best'.
%                   [SaveOpts]      : save options to a global variable called nr_search_opts. Default: false;
%                   [ReuseOpts]     : reuse global options nr_search_opts. Default: false;
%
% OUT:
%    count     : 1*L vector of number of neighbors
%    idxs      : L cell array (or L times NN array if opts.ReturnAsArrays==1) of neighbor indices
%    dists     : L cell array (or L times NN array if opts.ReturnAsArrays==1) of distances of the neighbors
%
% (c) Copyright Duke University, 2008
% Mauro Maggioni
%
% EXAMPLE:
%   X=[cos(linspace(0,2*pi-2/500*pi,500));sin(linspace(0,2*pi-2/500*pi,500))];
%   [count,idxs, dists,NNInfo] = nrsearch(X, X, 3, 0)
%   [count,idxs, dists,NNInfo] = nrsearch(X, X, 3, 0,[],struct('ReturnAsArrays',0));
%   [count,idxs, dists,NNInfo] = nrsearch(X, (1:3)', 3, 0,[],struct('ReturnAsArrays',0));


count   = [];
idxs    = [];
dists   = [];
NNInfo  = [];

% Setup parameters
if nargin<6,    opts = []; end;
if nargin<5,    metric = []; end;
if (nargin<4) || isempty(radius),    radius = 0; end;

if isfield(opts,'ReuseOpts') && (opts.ReuseOpts),
    global nr_search_opts;
    nr_search_opts.ReuseOpts = opts.ReuseOpts;
    nr_search_opts.SaveOpts = opts.SaveOpts;
    opts = nr_search_opts;
else
    if ~isfield(opts,'DistRestrict'),           opts.DistRestrict = {}; end;
    if ~isfield(opts,'Tolerance'),              opts.Tolerance = 0; end;
    if ~isfield(opts,'ReturnDistSquared'),      opts.ReturnDistSquared = 0; end;
    if ~isfield(opts,'ReturnAsArrays'),         opts.ReturnAsArrays = 0; end;
    if ~isfield(opts,'DistInfo'),               opts.DistInfo = []; end;
    if ~isfield(opts,'NNInfo'),                 opts.NNInfo = []; end;
    if ~isfield(opts.NNInfo,'Atria'),           opts.NNInfo.Atria = []; end;
    if ~isfield(opts.NNInfo,'CoverTree'),       opts.NNInfo.CoverTree = []; end;
    if (~isfield(opts,'XIsTransposed')) || (isempty(opts.XIsTransposed)),       opts.XIsTransposed = false; end;
    if isempty(NN),    NN = 0; end;
    if (~isfield(opts,'FastNNSearcher')) || (isempty(opts.FastNNSearcher)),     opts.FastNNSearcher = 'best'; end;
    if ~isfield(opts,'SaveOpts'),               opts.SaveOpts = false; end;
    if ~isfield(opts,'ReuseOpts'),              opts.ReuseOpts = false; end;
end;

idxs = []; count = []; dists = []; NNInfo = [];

if isfield(opts,'DistInfo') && isempty(opts.DistInfo)
    if strcmpi(opts.FastNNSearcher,'best'),
        lUseCoverTrees  = (exist('covertree')==3);
        lUseTSTool      = (exist('nn_prepare')==3);
        lUseANNsearch   = (exist('ANNsearch')==3);
    elseif strcmpi(opts.FastNNSearcher,'nn_search'),
        lUseTSTool      = true;
        lUseANNsearch   = false;
        lUseCoverTrees  = false;
    elseif strcmpi(opts.FastNNSearcher,'ANNsearch'),
        lUseTSTool      = false;
        lUseANNsearch   = true;
        lUseCoverTrees  = false;
    elseif strcmpi(opts.FastNNSearcher,'covertree'),
        lUseTSTool      = false;
        lUseANNsearch   = false;
        lUseCoverTrees  = true;    
    else
        lUseTSTool      = false;
        lUseANNsearch   = false;
        lUseCoverTrees  = false;
    end;
else
    lUseTSTool      = false;
    lUseANNsearch   = false;
    lUseCoverTrees  = false;
end;

if (lUseTSTool==false) && (lUseANNsearch==false) && (lUseCoverTrees==false),
    warning('\n ***nrsearch: could not find any fast nearest neighbor package installed. See ''help nrsearch''. Aborting nearest neighbor search.');
    return;
end;

% Find dimension and number of points
if opts.XIsTransposed==false,
    [lDim,lNPts] = size(cX); lNQueryPts = size(cXq,2);
else
    [lNPts,lDim] = size(cX); lNQueryPts = size(cXq,1);
end;

% If no cXq, then query points are all points
if isempty(cXq),    cXq = uint32(1:lNPts);  end;
lGoWithIdxs = (strcmpi(class(cXq),'uint32'));

lUseCoverTrees = false;

% A fast nearest neighbor searcher is available
% Check if a set of candidate nearest neighbors is provided
if isempty(opts.DistRestrict),
    if lUseANNsearch || lUseCoverTrees,
        % Use ANNsearch
        if opts.XIsTransposed==true,
            cX = cX';
            cXq = cXq';
        end;
        lNumberOfPoints = size(cX,2);
        % If cXq is a set of indices into cX, unfortunately here we need to create it as a set of points
        if lGoWithIdxs, cXq = cX(:,cXq); end;
        if nargout>2,
            if radius==0,
                if ~lUseCoverTrees,
                    [idxs,dists] = ANNsearch(cX, cXq, NN, opts.Tolerance);
                    if ~opts.ReturnDistSquared,
                        for k = 1:size(dists,1);dists(k,:) = sqrt(dists(k,:));end;
                    end;
                else
                    if (isempty(opts.NNInfo.CoverTree))
                        A.vectors=cX;A.theta=.5;A.maxdescend=int32(round(log2(lNPts)));NNInfo.CoverTree=covertree(A);
                    else
                        NNInfo.CoverTree = opts.NNInfo.CoverTree;
                    end;
                    [idxs,dists] = covertree_nn_search(cX,NNInfo.CoverTree,cXq,NN);idxs=idxs';dists=dists';
                end;
                dists = dists';
                idxs  = idxs';
                if ~opts.ReturnAsArrays,
                    idxs   = mat2cell(idxs,ones(size(idxs,1),1),size(idxs,2));
                    dists  = mat2cell(dists,ones(size(dists,1),1),size(dists,2));
                end;
            else
                if ~lUseCoverTrees,
                    [count, idxs, dists] = ANNrsearch(cX, cXq, NN, radius^2, opts.Tolerance);
                    if ~opts.ReturnDistSquared,
                        for k = 1:length(dists);dists{k} = sqrt(dists{k});end;
                    end;
                else
                    if (isempty(opts.NNInfo.CoverTree))
                        A.vectors=cX;A.theta=.5;A.maxdescend=int32(round(log2(lNPts)));NNInfo.CoverTree=covertree(A);
                    else
                        NNInfo.CoverTree = opts.NNInfo.CoverTree;
                    end;
                    [idxs,dists] = covertree_range_search(cX,NNInfo.CoverTree,cXq,radius);
                end;
                if opts.ReturnAsArrays,
                    idxs_2(lNumberOfPoints,lNumberOfPoints) = uint32(0);
                    dists_2(lNumberOfPoints,lNumberOfPoints) = single(Inf);
                    for k = 1:length(idxs),
                        idxs_2(k,1:length(idxs{k})) = uint32(idxs{k});
                        dists_2(k,1:length(dists{k})) = dists{k};
                    end;
                    idxs = idxs_2; dists = dists_2; clear idxs_2 dists_2;
                end;
            end;
        else
            % Possibly save some computations and memory by not computing dists
            if radius==0,
                if ~lUseCoverTrees,
                    [idxs] = ANNsearch(cX, cXq, NN, opts.Tolerance);
                else
                    if (isempty(opts.NNInfo.CoverTree))
                        A.vectors=cX;A.theta=.5;A.maxdescend=int32(round(log2(lNPts)));NNInfo.CoverTree=covertree(A);
                    else
                        NNInfo.CoverTree = opts.NNInfo.CoverTree;
                    end;
                    [idxs] = covertree_nn_search(cX,NNInfo.CoverTree,cXq,NN);                    
                end;
                if ~opts.ReturnAsArrays,
                    idxs  = uint32(mat2cell(idxs,ones(size(idxs,1),1),size(idxs,2)));
                end;
            else
                if ~lUseCoverTrees,
                    [count, idxs] = ANNrsearch(cX, cXq, NN, radius^2, opts.Tolerance);
                    idxs = uint32(idxs);
                else
                    if (isempty(opts.NNInfo.CoverTree))
                        A.vectors=cX;A.theta=.5;A.maxdescend=int32(round(log2(lNPts)));NNInfo.CoverTree=covertree(A);
                    else
                        NNInfo.CoverTree = opts.NNInfo.CoverTree;
                    end;
                    [idxs,dists] = covertree_range_search(cX,NNInfo.CoverTree,cXq,radius);
                end
            end;
        end;
    elseif lUseTSTool,
        % Use nn_search
        if opts.XIsTransposed==false,
            lX = cX';
            if ~lGoWithIdxs, lXq = cXq'; else lXq=cXq; end;
        else
            lX = cX;
            lXq = cXq;
        end;
        % Unfortunately there's a bug in nn_search, when there is only one point
        if size(lX,1)==1,
            lOnlyOnePoint = true;
            lX = [lX;lX];
        else
            lOnlyOnePoint = false;
        end;
        if (isempty(opts.NNInfo.Atria)) || (opts.NNInfo.Atria.params(1)~=size(lX,1)),
            NNInfo.Atria = nn_prepare( lX );
        else
            NNInfo.Atria = opts.NNInfo.Atria;
        end;
        lNumberOfPoints = size(lX,1);
        if radius==0,
            % Perform nearest neighbors search
            NN = min(size(lX,1),NN);
            [idxs,dists] = nn_search( lX, NNInfo.Atria, double(lXq), NN,-1,opts.Tolerance );
            idxs = uint32(idxs);
            if lOnlyOnePoint,
                idxs = ones(size(idxs));
            end;
            % Format output as requested
            if opts.ReturnDistSquared,
                for k = 1:size(dists,1);dists(k,:) = dists(k).^2;end;
            end;
            if ~opts.ReturnAsArrays,
                idxs  = mat2cell(idxs,ones(size(idxs,1),1),size(idxs,2));
                dists  = mat2cell(dists,ones(size(dists,1),1),size(dists,2));
            end;
            count = NN;
        else
            % Perform range search
            [count, neighbors] = range_search( lX, NNInfo.Atria, double(lXq), radius,-1 );
            % The results are NOT sorted. Sort them
            for p = 1:size(neighbors,1),
                [neighbors{p,2},lSortedIdxs] = sort(neighbors{p,2});
                neighbors{p,1} = neighbors{p,1}(lSortedIdxs);
            end;
            % Convert to standard output format
            if ~lOnlyOnePoint,
                if NN==0,
                    for k = length(count):-1:1,
                        idxs{k} = uint32(neighbors{k,1});
                        dists{k} = neighbors{k,2};
                    end;
                else
                    % Truncate to nearest NN points
                    for k = length(count):-1:1,
                        count(k) = min(count(k),NN);
                        idxs{k} = uint32(neighbors{k,1}(1:count(k)));
                        dists{k} = neighbors{k,2}(1:count(k));
                    end;
                end;
            else
                for k = length(count):-1:1,
                    idxs{k} = uint32(1);
                    dists{k} = neighbors{k,2};
                end;
            end;
            clear neighbors;
            if opts.ReturnDistSquared,
                for k = 1:length(count);dists{k} = dists{k}.^2;end;
            end;
            if opts.ReturnAsArrays,
                idxs_2 (lNumberOfPoints,lNumberOfPoints) = uint32(0);
                dists_2(lNumberOfPoints,lNumberOfPoints) = single(Inf);
                for k = 1:length(idxs),
                    idxs_2(k,1:length(idxs{k})) = uint32(idxs{k});
                    dists_2(k,1:length(dists{k})) = dists{k};
                end;
                idxs = uint32(idxs_2); dists = dists_2; clear idxs_2 dists_2;
            end;
        end;
    end;
else
    % A list of candidate nearest neighbors was provided
    count = zeros(lNPts,1);
    for k = lNPts:-1:1,
        if (iscell(opts.DistRestrict)) && (length(opts.DistRestrict)==lNPts),
            % A NN list is provided for each point
            lCurIdxs = opts.DistRestrict{k};
        else
            % Only one NN list, interpret as relative to the current point
            lCurIdxs = k+opts.DistRestrict;
            % Clip indices out of range
            lCurIdxs(find((lCurIdxs<1) | (lCurIdxs>lNPts))) = [];
        end;
        % Find the neighbors
        if radius==0,
            if ~lUseTSTool,
                if opts.XIsTransposed,
                    cX = cX';
                end;
                [idxs_tmp,dists_tmp] = ANNsearch( cX(:,lCurIdxs), cX(:,k), NN, opts.Tolerance );
                if ~opts.ReturnDistSquared,
                    % Compute the sqrt of the distances
                    dists_tmp = sqrt(dists_tmp);
                end;
            else
                if opts.XIsTransposed==false,
                    lX = cX(:,lCurIdxs)';
                else
                    lX = cX(lCurIdxs,:);
                end;
                if isempty(opts.NNInfo.Atria),          % MM: Not correct if Atira is not for the restricted indices
                    NNInfo.Atria = nn_prepare( lX );
                else
                    NNInfo.Atria = opts.NNInfo.Atria;
                end;
                [idxs_tmp,dists_tmp] = nn_search( lX, NNInfo.Atria, cX(:,k), NN, opts.Tolerance );
                if opts.ReturnDistSquared,
                    for k = 1:size(idxs_tmp,1);dists_tmp(k,:) = dists_tmp(k,:).^2;end;
                end;
            end;
        else
            if ~lUseTSTool,
                if opts.XIsTransposed,
                    cX = cX';
                end;
                [count_tmp,idxs_tmp,dists_tmp] = ANNrsearch( cX(:,lCurIdxs), cX(:,k), NN, radius^2, opts.Tolerance );
                if ~opts.ReturnDistSquared,
                    % Compute the sqrt of the distances
                    dists_tmp = sqrt(dists_tmp);
                end;
            else
                if opts.XIsTransposed==false,
                    lX = cX(:,lCurIdxs)';
                else
                    lX = cX(lCurIdxs,:);
                end;
                if isempty(opts.NNInfo.Atria) || (opts.NNInfo.Atria.params(1)~=size(lX,1)),
                    NNInfo.Atria = nn_prepare( lX );
                else
                    NNInfo.Atria = opts.NNInfo.Atria;
                end;
                [count_tmp, neighbors_tmp] = range_search( lX, NNInfo.Atria, cX(k,:)', radius );
                for k = length(count_tmp):-1:1,
                    idxs_tmp{k} = neighbors_tmp{k,1};
                    dists_tmp{k} = neighbors_tmp{k,2};
                end;
                clear neighbors;
            end;
        end;
    end;
    if ~opts.ReturnAsArrays,
        idxs{k} = uint32(lCurIdxs(idxs_tmp));
        dists{k} = dists_tmp;
    else
        idxs(k,:) = uint32(lCurIdxs(idxs_tmp));
        if argout>2,
            dists(k,:) = dists_tmp;
        end;
    end;
end;

% Build count
if (~iscell(idxs)) && (size(idxs,1)==1),
    count = ones(1,length(idxs),'uint32');
else
    if ~opts.ReturnAsArrays,
        for k = length(idxs):-1:1,
            count(k) = length(idxs{k});
        end
    else
        count = size(idxs,2)*ones(1,size(idxs,1));
    end;
end;

if opts.SaveOpts,
    global nr_search_opts
    if exist('NNInfo'),
        opts.NNInfo = NNInfo;
    end;
    nr_search_opts = opts;
end;


return;