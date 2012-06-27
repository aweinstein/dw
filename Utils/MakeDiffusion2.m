function [T, W, P, DistInfo] = MakeDiffusion2(Points, Radius, Options)

%
% function [T, W, P] = MakeDiffusion2(Points, Radius, Options)
%
% MAKEDIFFUSION2 builds a diffusion operator on a cloud of points embedded in
% a metric space.  A matrix of graph weights for the point cloud is
% constructed by building a local bump function centered at each point and
% then the weight matrix is normalized to form the diffusion operator T.
%
% MAKEDIFFUSION2 builds the weight matrix by constructing each column as a local
% bump function.  For each point x, the points within Radius arre found. The
% weight
%
%    exp( - d(x,y) * d(y,x) / (epsilon_x * epsilon_y) )
%
% is then assigned to the pair (x,y).  Note that because of autotuning,
% d(x,y) is not necessarily equal to d(y,x) ... thus the product.  Also
% note that the constant epsilon can vary from point to point.
%
% The matrix of weights is then normalized, in one of a number of possible
% ways, to form a final "diffusion" operator T.  Normally, either a symmetric
% Markov or symmetric Laplace-Beltrami normalization is used in order to make
% eigensolving
%
%
% In:
%    Points    = DxN matrix specifying N points in R^D
%    Radius    = radius of range search around each point. Can be a vector. If Radius==0, uses kNN only.
%    Options   = structure for specifying algorithm options which can contain
%                zero or more of the following fields:
%
%      Display        : when this is zero the only output will be error messages
%      Epsilon        : either a scalar or a vector of length N giving time-step constants, to construct N corresponding different graphs.
%      [kNN]              : truncate the neighbors to the k closest points.
%      [kNNAutotune]      : rescale Epsilon based on the distance to the kNNAutotune nearest point
%      Normalization  : method for normalizing the matrix of weights
%         'bimarkov'        force row and column sums to be 1
%         'beltrami'        Laplace-Beltrami normalization ala Coifman-Lafon
%         'markov'          force row sums to be 1
%         'sbeltrami'       symmetric conjugate to beltrami
%         'smarkov'         symmetric conjugate to markov
%         'FokkerPlanck'    Fokker-Planck normalization
%         'sFokkerPlanck'   symmetric conjugate to Fokker-Planck normalization
%      [Distance]      : distance to use
%         'Euclidean'   : default
%         A combination of the following strings (case sensitive):
%         'PCA'                 : project onto local tangent plane
%         'direct' or 'inverse' : local pca or its inverse
%         'cov'                 : normalize by number of local points
%         'noS'                 : no scaling by the singular values
%      [PCARadius]     : radius for the computation of local PCA. it can be a vector (matching Radius). Default: Radius.
%      [PCAtruncD]     : only use the first D principal axes for computing local PCA distance. Default: D is maximum possible.
%
%    [Metric]           : handle to a metric function
%    [DistInfo]         : if not empty, it will be used for distance computation. It is assumed that it contains at least
%                         the neighbors up to the requested distance Radius (and PCARadius if provided), or kNN (if provided).
%                         It assumes the neighbors are listed in increasing distance.
%    [Symmetrization]   : how to symmetrize the weight matrix:
%                         'W.Wt'     : replaces W by W.*W'. This is the default value
%                         'WWt'      : replaces W by W*W'
%                         'WtW'      : replaces W by W'*W
%                         'W+Wt'     : replaces W by W+Wt
%
%    The default options are:
%
%       Options.Display = 1;
%       Options.Epsilon = 1.0;
%       Options.Normalization = 'beltrami';
%
% Out:
%    T         = NxN sparse matrix giving the normalized diffusion operator
%    W         = NxN sparse matrix of weights
%    P         = in the case of symmetric normalizations, this is the NxN diagonal
%                matrix which relates the nonsymmetric version to the symmetric
%                form via conjugation; symmetric = D^(-1)*nonsymetric*D
%    DistInfo  = structure containing distance info:
%                   count    :   N vector of counts of points near each point
%                   idxs     :   N cell of indices of the nearest points to each point
%                   dists    :   N cell of distances between the nearst points and each point
%
%
% Example:
%   X=[cos(linspace(0,2*pi-0.0001,500));sin(linspace(0,2*pi-0.0001,500))];
%   T= MakeDiffusion2(X, 0.05, struct('Normalization','bimarkov'));

% Author:
%     jcb   8/06 last revision
%     mm    4/07 some small changes and bug fixes
%     mm    9/9  added F-P normalization option
%     mm    10/18 small changes and bug fixes
%     mm    11/18
%     mm    12/10
%
%
% (c) Copyright Yale University, 2006, James C Bremer Jr.
% (c) Copyright Duke University, 2007, Mauro Maggioni
%

N = size(Points,2);

% initialize output variables in case of error return
T = [];
W = [];
D = [];

%
% process options
%
Display = 1;
Epsilon = 1.0;
Normalization = 'beltrami';
Distance = 'Euclidean';
PCARadius = Radius;
PCAtruncD = Inf;
DistInfo = [];
kNN = Inf;
kNNAutotune = 0;
Symmetrization = 'W.Wt';

lDistInfoProvided = false;

if nargin<3,
    Options = struct();
end

fn = fieldnames(Options);
for j=1:length(fn)
    name = fn{j};
    value = getfield(Options, name);

    if strcmpi(name, 'Display')
        Display = value;
    elseif strcmpi(name, 'Epsilon')
        Epsilon = value;
    elseif strcmpi(name, 'Normalization')
        Normalization = value;
    elseif strcmpi(name, 'Metric')
        Metric = value;
    elseif strcmpi(name,'Distance'),
        Distance = value;
    elseif strcmpi(name,'PCARadius'),
        PCARadius = value;
    elseif strcmpi(name,'PCAtruncD'),
        PCAtruncD = value;
    elseif strcmpi(name,'DistInfo'),
        lDistInfoProvided = true;
        DistInfo = value;
    elseif strcmpi(name,'kNN'),
        kNN = value;
    elseif strcmpi(name,'kNNAutotune'),
        kNNAutotune = value;
    elseif strcmpi(name,'Symmetrization'),
        Symmetrization = value;
    else
        fprintf('MakeDiffusion.m: invalid option "%s" ignored.\n', name);
    end
end

if Display > 1
    fprintf('Options: \n');
    fprintf('\tDisplay = %d\n', Display);
    fprintf('\tEpsilon = %g\n', Epsilon);
    fprintf('\tNormalization = %s\n', Normalization);
end

%
% Perform range/nearest neighbor searches on the input set
%
if Display; fprintf('performing range searches ... '); TIME = clock; end;

if isempty(DistInfo),
    % If Radius>0, do a range search
    if Radius>0,
        if ~exist('Metric')
            [DistInfo.count, DistInfo.idxs, DistInfo.dists] = rsearch(Points, [], max([Radius,PCARadius]));
        else
            [DistInfo.count, DistInfo.idxs, DistInfo.dists] = rsearch(Points, [], max([Radius,PCARadius]), Metric);
        end;    
    else % otherwise use nearest neighbors
        if ~exist('Metric')
            [DistInfo.idxs, DistInfo.dists] = nnsearch(Points, [], kNN);
        else
            [DistInfo.idxs, DistInfo.dists] = nnsearch(Points, [], kNN, Metric);
        end;
        DistInfo.count = kNN;
    end;
    save DistInfo DistInfo
end;

if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;

N = size(Points,2);

%
% Compute the weights and associated operators
% 
if Display; fprintf('computing weights...'); TIME = clock; end;
for k = 1:length(Radius),
    % Find the neighbors of each point
    for i = N:-1:1,
        if (Radius>0) & ((length(Radius)>1) | (lDistInfoProvided==true)),
            % If a Radius>0 is specified, and either there's more than one radius or DistInfo was provided,
            % restrict the neighbors to Radius(k) and the specified kNN.
            try
                idxs{i} = find(DistInfo.dists{i}<=Radius(k));
            catch
                fprintf('Ops!');
            end;
            if length(idxs{i})>kNN,
                idxs{i} = idxs{i}(1:kNN);
            end;
            dists{i} = DistInfo.dists{i}(idxs{i});
            idxs{i}  = DistInfo.idxs{i}(idxs{i});
            count(i) = length(idxs{i});
        else
            % If Radius==0, or if there is only one Radius and DistInfo was not provided
            dists{i} = DistInfo.dists{i};
            idxs{i}  = DistInfo.idxs{i};
            count(i) = length(idxs{i});
        end;
    end;
    % Do the autotuning of distances if requested
    if kNNAutotune > 0
        if Display; fprintf('autotuning distances ... '); TIME = clock; end;
        for j=1:N
            temp = sort(dists{j}, 'ascend');              % MM: This could be once, at the cost of a bit of memory.
            lMaxTempIdxs = min([kNNAutotune,length(temp)]);
            if temp(lMaxTempIdxs)==0,
                dists{j} = 0;
            else
                dists{j} = dists{j}./temp(lMaxTempIdxs);
            end;
        end
        if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;
    end
    % Do the local PCA warping of the distances
    if ~isempty(strfind(Distance,'PCA')),
        % Compute the local PCA matrices
        for i = N:-1:1,
            % Find the nearest neighbors to be used for the local PCA computation
            if PCARadius(k)<max([Radius,PCARadius]),
                pca_idxs  = DistInfo.idxs{i}(find(DistInfo.dists{i}<=PCARadius(k)));
                pca_count = length(pca_idxs{i});
            else
                pca_idxs  = DistInfo.idxs{i};
                if iscell(DistInfo.count),
                    pca_count = DistInfo.count{i};
                else
                    pca_count = DistInfo.count(i);
                end;
            end;
            % Center the points
            lCenteredPcaPts = Points(:,pca_idxs)-repmat(Points(:,i),[1,length(pca_idxs)]);
            % Compute local PCA
            [lU{i},lS{i}] = svd(lCenteredPcaPts,0);lS{i}=diag(lS{i});
            % Fix case in which there were less points than the dimension
            if size(lCenteredPcaPts,2)<size(lCenteredPcaPts,1),         % If there are less points than dimensions
                lU{i} = [lU{i},zeros(size(lU{i},1),size(lCenteredPcaPts,1)-size(lCenteredPcaPts,2))];
                lS{i} = [lS{i};zeros(size(lCenteredPcaPts,1)-size(lCenteredPcaPts,2),1)];
            end;
            % Truncate the PCAs if requested
            if PCAtruncD<length(lS{i}),
                lS{i}(PCAtruncD+1:length(lS{i}))=0;
            end;
            % Compute distances between points based on the local PCA
            lCenteredRadiusPts = Points(:,idxs{i})-repmat(Points(:,i),[1,length(idxs{i})]);
            if ~isempty(strfind(Distance,'cov')),
                lS{i} = lS{i}/double(pca_count);
            end;
            if ~isempty(strfind(Distance,'noS')),
                dists{i} = sqrt(sum((lU{i}*lCenteredRadiusPts).^2,1)');
            elseif ~isempty(strfind(Distance,'direct')),
                dists{i} = sqrt(sum((diag(lS{i})*lU{i}*lCenteredRadiusPts).^2,1)');
            elseif ~isempty(strfind(Distance,'inverse')),
                try
                    dists{i} = sqrt(sum((diag(lS{i}.^(-1))*lU{i}*lCenteredRadiusPts).^2,1)');
                catch
                    fprintf('Something wrong.');
                end;
            end;
        end;
    end;
    if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;

    % Build the weight matrix
    if Display; fprintf('forming weight matrix ... '); TIME = clock; end;
    % First form the distance matrix
    rowidxs = zeros(sum(count), 1);
    colidxs = zeros(sum(count), 1);
    distsmat = zeros(sum(count), 1);
    index = 1;
    location = 1;
    for j=1:length(count),
        rowidxs(location:(location+count(j)-1))     = idxs{index};
        colidxs(location:(location+count(j)-1))     = index;
        distsmat(location:(location+count(j)-1))    = dists{index};
        location                                    = location+count(j);
        index                                       = index+1;
    end;
    
    if Epsilon(k)>0,
        W = sparse(double(rowidxs), double(colidxs),  distsmat, N, N, sum(count));
    else
        W = sparse(double(rowidxs), double(colidxs),  ones(size(distsmat)), N, N, sum(count));
    end;
    clear idxs colidxs dists count;

    % now take d(x,y) = d(x,y)*d(y,x) in order to symmetrize
    if strcmpi(Symmetrization,'W.Wt'),
        W = W .* W';
    elseif strcmpi(Symmetrization,'WWt'),
        W = W*W';
    elseif strcmpi(Symmetrization,'WtW'),
        W = W'*W;
    elseif strcmpi(Symmetrization,'W+Wt'),
        W = W + W';
    end;

    % adjust with epsilon
    if Epsilon(k)>0,
        D = sparse(1:N, 1:N, 1./Epsilon(k), N, N, N);
        W = D*W*D;
        % take exp(-W) ... be careful to add the diagonal back in
        [i j s] = find(W);
        nnzW=nnz(W);
        clear W;

        i = [i; (1:N)'];
        j = [j; (1:N)'];
        s = [s; zeros(N,1)];
        W = sparse(i, j, exp(-s), N, N, nnzW+N);
    else
        W = W + speye(size(W));
    end;

    if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;

    P = [];

    %
    % Form the diffusion operator
    %
    if ~strcmpi(Normalization, 'bimarkov')
        if Display; fprintf('normalizing to form operator '); TIME = clock; end;
    end

    if strcmpi(Normalization, 'bimarkov')
        T = Bimarkov(W);
    elseif strcmpi(Normalization, 'smarkov')
        if Display; fprintf('(symmetric markov) ... '); end;

        D = sparse(1:N, 1:N, 1./sqrt(sum(W,2)), N, N, N);
        P = D;
        T = D*W*D;

        T = (T'+T)/2;  % iron out numerical wrinkles

    elseif strcmpi(Normalization, 'markov')
        if Display; fprintf('(markov) ... '); end;
        T = sparse(1:N, 1:N, 1./sum(W,2), N, N, N) * W;

    elseif strcmpi(Normalization, 'sbeltrami')

        if Display;fprintf('(symmetric beltrami) ... '); end;

        % beltrami normalization
        P = sparse(1:N, 1:N, 1./sum(W,2), N, N, N);
        K = P*W*P;

        D = sparse(1:N, 1:N, 1./sqrt(sum(K,2)), N, N, N);
        P = D;
        T = D*K*D;

        T = (T'+T)/2;  % iron out numerical wrinkles


    elseif strcmpi(Normalization, 'beltrami')
        if Display; fprintf('(beltrami) ... '); end;

        % beltrami normalization
        D = sparse(1:N, 1:N, 1./sum(W,2), N, N, N);
        K = D*W*D;

        V = sparse(1:N, 1:N, 1./sum(K,2), N, N, N);
        T = V*K;

    elseif strcmpi(Normalization, 'FokkerPlanck')
        if Display; fprintf('(FokkerPlanck) ... '); end;

        % beltrami normalization
        D = sparse(1:N, 1:N, sqrt(1./sum(W,2)), N, N, N);
        K = D*W*D;

        D = sparse(1:N, 1:N, 1./sum(K,2), N, N, N);
        T = D*K;

    elseif strcmpi(Normalization, 'sFokkerPlanck')
        if Display; fprintf('(sFokkerPlanck) ... '); end;

        % Fokker Planck normalization
        D = sparse(1:N, 1:N, sqrt(1./sum(W,2)), N, N, N);
        K = D*W*D;

        D = sparse(1:N, 1:N, sqrt(1./sum(K,2)), N, N, N);
        P = D;
        T = D*K*D;

        T = (T+T')/2;           % Iron numerically wrinkles

    else
        fprintf('\n\nMakeDiffusion2.m: unknown normalization\n');
        return;
    end

    if ~strcmpi(Normalization, 'bimarkov')
        if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;
    end;

    if length(Radius)>1,
        vT{k} = T;
        if nargout>1,
            vW{k} = W;
            if nargout>2,
                vP{k} = P;
            end;
        end;
    end;

end;

if length(Radius)>1,
    T=vT;
    if nargout>1,
        W = vW;
        if nargout>2,
            P = vP;
        end;
    end;
end;

return;