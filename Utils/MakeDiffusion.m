function [T, W, P] = MakeDiffusion(Points, NN, Options)

%
% function [T, W, P] = MakeDiffusion(Points, NN, Options)
%
% MAKEDIFFUSION builds a diffusion operator on a cloud of points embedded in
% a metric space.  A matrix of graph weights for the point cloud is
% constructed by building a local bump function centered at each point and
% then the weight matrix is normalized to form the diffusion operator T.
%
% MAKEDIFFUSION builds the weight matrix by constructing each column as a local
% bump function.  For each point x, the closest NN neighbors are found.  If
% the autotune option is specified, the distances are adjusted so that the
% Options.Autotune^th neighbor is at a distance 1.0 from the point x.  The
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
%    NN        = number of nearest neighbors
%    Options   = structure for specifying algorithm options which can contain
%                zero or more of the following fields:
%
%      Display        : when this is zero the only output will be error
%                       messages
%      Autotune       : controls autotuning, as above. Can be a vector.
%      Epsilon        : either a scalar or a vector of length N giving
%                       time-step constants
%      Normalization  : method for normalizing the matrix of weights
%         'bimarkov'    force row and column sums to be 1
%         'beltrami'    Laplace-beltrami normalization ala Coifman-Lafon
%         'markov'      force row sums to be 1
%         'sbeltrami'   symmetric conjugate to beltrami
%         'smarkov'     symmetric conjugate to markov
%         'FokkerPlanck'Fokker-Planck normalization
%
%    Metric           : handle to a metric function
%    Mod              : when computing distances, mod them in [0,Mod]
%
%    The default options are:
%
%       Options.Display = 1;
%       Options.Autotune = 3;
%       Options.Epsilon = 1.0;
%       Options.Normalization = 'beltrami';
%
% Out:
%    T         = NxN sparse matrix giving the normalized diffusion operator
%    W         = NxN sparse matrix of weights
%    P         = in the case of symmetric normalizations, this is the NxN diagonal
%                matrix which relates the nonsymmetric version to the symmetric
%                form via conjugation; symmetric = D^(-1)*nonsymetric*D
%
%
% Example:
%   X=[cos(linspace(0,2*pi-0.0001,500));sin(linspace(0,2*pi-0.0001,500))];
%   T= MakeDiffusion(X, 3, struct('Normalization','bimarkov'));
%
% Author:
%     jcb   8/06 last revision
%     mm    4/07 some small changes and bug fixes
%     mm    9/9  added F-P normalization option
%
% (c) Copyright Yale University, 2006, James C Bremer Jr.
% (c) Copyright Duke University, 2007, Mauro Maggioni
%

N = size(Points,2);

NN = min(NN, N);

% initialize output variables in case of error return
T = [];
W = [];
D = [];

%{
% process options
%}
Display = 1;
Epsilon = 1.0;
Autotune = 3;
Normalization = 'beltrami';

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
    elseif strcmpi(name, 'Autotune')
        Autotune = value;
    elseif strcmpi(name, 'Normalization')
        Normalization = value;
    elseif strcmpi(name, 'Metric')
        Metric = value;
    else
        fprintf('MakeDiffusion.m: invalid option "%s" ignored.\n', name);
    end
end

if Display > 1
    fprintf('Options: \n');
    fprintf('\tNN = %d\n', NN);
    fprintf('\tDisplay = %d\n', Display);
    %fprintf('\tEpsilon = %g\n', Epsilon);
    fprintf('\tAutotune = %d\n', Autotune);
    fprintf('\tNormalization = %s\n', Normalization);
end

% adjust autotune if necessary
Autotune = min(Autotune, NN-1);

% if necessary, make Epsilon a vector for future convienence
if length(Epsilon)==1
    Epsilon = repmat(Epsilon, N, 1);
end

%{
% perform nearest neighbors searches on the input set
%}
if Display; fprintf('performing nn searches ... '); TIME = clock; end;

if ~exist('Metric')
    [idxs, dists] = nnsearch(Points, Points, NN,[],struct('ReturnAsArrays',true));
    idxs=idxs';dists=dists';
else
    [idxs, dists] = nnsearch(Points, Points, NN, Metric,struct('ReturnAsArrays',true));
end

if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;


for lk = 1:length(Autotune),
    % find the autotuning constants
    if Autotune(lk) > 0
        cautotune = zeros(N, 1);      % autotune constants
        if Display; fprintf('autotuning distances ... '); TIME = clock; end;
        for j=1:N
            temp = sort(dists(:,j), 'ascend');
            if temp(Autotune(lk)+1)==0,
                dists(:,j) = 0;
            else
                dists(:,j) = dists(:,j)./temp(Autotune(lk)+1);
            end;
        end
        if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;
    end

    %{
    % form the weights matrix in stages
    %}
    if Display; fprintf('forming weight matrix ... '); TIME = clock; end;

    % first form the distance matrix
    %W = sparse(idxs, diag(1:N)*ones(size(idxs)), dists, N, N, NN*N);
    %W = sparse(double(idxs), kron((1:N)', ones(1,NN)), dists, N, N, NN*N);

    W = sparse(double(idxs), kron(ones(NN,1), (1:N)), dists, N, N, NN*N);

    % now take d(x,y) = d(x,y)*d(y,x) in order to symmetrize
    W = W .* W';

    % adjust with epsilon
    D = sparse(1:N, 1:N, 1./Epsilon, N, N, N);
    W = D*W*D;

    % take exp(-W) ... be careful to add the diagonal back in
    [i j s] = find(W);

    i = [i; (1:N)'];
    j = [j; (1:N)'];
    s = [s; zeros(N,1)];
    W = sparse(i, j, exp(-s), N, N, nnz(W)+N);
    if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;

    P = [];

    %{
    % form the operator
    %}
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

        V = sparse(1:N, 1:N, 1./sum(K,2), N, N, N);
        T = V*K;

    elseif strcmpi(Normalization, 'sFokkerPlanck')
        if Display; fprintf('(sFokkerPlanck) ... '); end;

        % beltrami normalization
        D = sparse(1:N, 1:N, sqrt(1./sum(W,2)), N, N, N);
        K = D*W*D;

        V = sparse(1:N, 1:N, sqrt(1./sum(K,2)), N, N, N);
        T = V*K*V;

        T = (T+T')/2;           % Iron numerically wrinkles

    else
        fprintf('\n\nMakeDiffusion2.m: unknown normalization\n');
        return;
    end

    if ~strcmpi(Normalization, 'bimarkov')
        if Display; fprintf('%g seconds\n', etime(clock, TIME)); end;
    end
end;

return;