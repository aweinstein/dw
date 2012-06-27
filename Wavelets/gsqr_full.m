function [Idxs,DIdxs,Q,R,RI] = gsqr_full(A, Options)

% function [Idxs DIdxs Q R RI] = gsqr_full(A, [Options])
%
% GSQR computes a partial QR decomposition of the matrix A via modified
% Gram-Schmidt orthogonalization.  A rank revealing QR decomposition is a
% factorization
%
%      A*P = | Q11 Q12 | * | R11 R12 |
%                          |  0  R22 |
%
% where P is a permutation, Q = [Q11 Q12] is orthogonal, and R22 is a matrix
% with small L^2 norm.  Under these assumptions, the contribution to the
% matrix product from the term Q12*R22 can be safely ignored (introducing
% a controllable error).
%
% GSQR computes the matrices Q11 and R11 in the RRQR decomposition, which
% yields the approximate factorization
%
%      [A(:,Idxs) A(:,DIdxs)] = [Q11*R11 Q11*R12]
%
% where R12 = Q11'*A(:,DIdxs).
%
% GSQR can also compute the inverse, IR11, of the matrix R11, which can be
% used to form a matrix
%
%      Proj = IR11*R12
%
% such that the diagram
%                  A
%        R^N ------------> R^M
%         |                / \
%    Proj |                 |
%         |                 |  T(:,Idxs)
%        \ /                |
%        R^K ------------> R^K
%              Identity
%
% is approximately commutative.
%
% If the returns variables R or RI are not specified then they will not be
% computed.  This can save a significant about of time.
%
% This is a matlab reference implementation of GSQR; its behavior is (almost)
% identical to the C MEX version but it isn't required to be efficient, only
% easy to understand and maintain.
%
% In:
%   A        = Input matrix
%   Options  = A structure specifying algorithm settings.
%
%   (1) Stopping conditions
%
%     StopPrecision   : stop when the l^2 norm of all remaining columns falls
%                       below this value
%     StopN           : specifies the maximum number of iterations
%     StopDensity     : stop when the density of all remaining columns becomes
%                       larger than this
%
%   (2) Thresholding
%
%     IPThreshold     : threshold for an inner product to be considered nonzero
%     Threshold       : threshold for the matrix Q
%
%   (3) Other Options
%
%     Reorthogonalize : perform two projections instead of one
%     Quiet           : when true, all output other than error messages is
%                       suppresed
%     ForceSparse     : force the working matrices Q, R, and RI to be
%                       sparse
%
%   The default options are as follows:
%
%     Options.StopPrecision      = eps;
%     Options.StopN              = Inf;
%     Options.StopDensity        = 1.0;
%     Options.IPThreshold        = 0.0;
%     Options.Reorthogonalize    = true;
%     Options.Quiet              = false;
%
% Out:
%   Idxs     = List of chosen columns, in the order in which they were chosen
%   DIdxs    = List of discarded columns, in arbitrary order
%   Q        = The matrix Q11 in the RRQR decomposition
%   R        = The matrix R11 in the RRQR decomposition
%   RI       = An approximate inverse for the matrix R
%
% Dependences:
%   none
%
% Version History:
%   jcb        10/2005        initial version based on FastGSOGP
%

T = cputime;

if ~exist('Options')

   Options = [];
end

if ~isfield(Options, 'StopN')
   Options.StopN = Inf;
end
if ~isfield(Options, 'StopPrecision')
   Options.StopPrecision = eps;
end
if ~isfield(Options, 'StopDensity')
   Options.StopDensity = 1.0;
end
if ~isfield(Options, 'Reorthogonalize')
   Options.Reorthogonalize = true;
end
if ~isfield(Options, 'IPThreshold')
   Options.IPThreshold = 0.0;
end
if ~isfield(Options, 'Threshold')
   Options.Threshold = 0.0;
end
if ~isfield(Options, 'ComputeR')
   Options.ComputeR = true;
end
if ~isfield(Options, 'ComputeRI')
   Options.ComputeRI = false;
end
if ~isfield(Options, 'Quiet')
   Options.Quiet = false;
end

if ~isfield(Options, 'ForceSparse')
    Options.ForceSparse = false;
end

% shortcut variables for options
StopN              = Options.StopN;
StopPrecision      = Options.StopPrecision;
StopDensity        = Options.StopDensity;
Reorthogonalize    = Options.Reorthogonalize;
IPThreshold        = Options.IPThreshold;
Threshold          = Options.Threshold;
Quiet              = Options.Quiet;
ForceSparse        = Options.ForceSparse;

ComputeR = nargout > 3;
ComputeRI = nargout > 4;


% setup some basic options
[M N] = size(A);

if ForceSparse
    Q = sparse([], [], [], M, N, N);
else
    Q = zeros(N, N);
end

if ComputeR
    if ForceSparse
        R = sparse([], [], [], N, N, N);
    else
        R = zeros(N,N);
    end
else
    R = [];
end

if ComputeRI
    if ForceSparse
        RI = sparse([], [], [], N, N, N);
    else
        RI = zeros(N,N);
    end
else
    RI = [];
end

Chosen = zeros(1,N);
Norms = zeros(1,N);
Idxs = zeros(1,N);
NumChosen = 0;

if Reorthogonalize
    NumProjections = 2;
else
    NumProjections = 1;
end

% compute norms
for j=1:N
    Norms(j) = norm(A(:,j));

    if Norms(j) == 0.0
      fprintf('norm 0\n');
     end

end

if ~Quiet
   fprintf('00000');
end

NoIps = 0;

for j=1:min(StopN, N);

    % choose a column
    [HeapNorm ChosenColumn] = max(Norms);

    Q(:,j) = A(:,ChosenColumn);
    A(:,ChosenColumn) = zeros(size(A,1), 1);

    ComputedNorm = norm(Q(:,j));

    for i=1:NumProjections

        % compute and threshold inner products
        ips = Q(:,1:j-1)'*Q(:,j);
        ips = ips .* (abs(ips) > IPThreshold*ComputedNorm);

        Q(:,j) = Q(:,j) - Q(:,1:j-1)*ips;

        if ComputeR
            R(1:j-1,j) = R(1:j-1,j)+ips;
        end
        if ComputeRI
            RI(:,j) = RI(:,j) - RI(:,1:j-1)*ips;
        end
    end

    % normalize Q
    ChosenNorm = norm(Q(:,j));
    Q(:,j) = Q(:,j)/ChosenNorm;

    if ComputeR
        R(j,j) = ChosenNorm;
    end

    if ComputeRI
        RI(j,j) = 1;
        RI(:,j) = RI(:,j)/ChosenNorm;
    end

    % threshold Q
    Q(:,j) = Q(:,j) .* (abs(Q(:,j)) > Threshold);

    if(ChosenNorm < StopPrecision)
      break;
    end

    if(nnz(Q(:,j))/M > StopDensity)
      break;
    end

    % mark column as chosen
    Chosen(ChosenColumn) = true;
    Norms(ChosenColumn) = -1.0;
    NumChosen = NumChosen+1;
    Idxs(NumChosen) = ChosenColumn;

    % update the norms of the overlapping columns
    ips = A'*Q(:,j);
    for i=find(abs(ips'./Norms)> IPThreshold)
        Norms(i) = sqrt(abs(Norms(i)^2 - ips(i)^2));
    end

    % display
    if ~Quiet && rand(1,1) < .25
      fprintf('\b\b\b\b\b%05d', j);
    end

% very verbose display
%    fprintf('column %d chosen with real norm %g, heap norm %g (difference %g) \n', ...
%      ChosenColumn, ChosenNorm, HeapNorm, HeapNorm-ChosenNorm);

end

% truncate list of indices and output matrix Q
Idxs = Idxs(1:NumChosen);
DIdxs = sort(find(~Chosen));
Q = sparse(Q(:, 1:NumChosen));
Q = Q .* (abs(Q) > Threshold);

% reorder R and RI
if ComputeR
   R = sparse(R(1:NumChosen, 1:NumChosen));
end

if ComputeRI
   RI = sparse(RI(1:NumChosen,1:NumChosen));
end

if ~Quiet
   if isfield(Options, 'dwtree')   % special gsqr mode
      fprintf('\b\b\b\b\b');
   else
      fprintf('\b\b\b\b\bdone (%d fcns chosen, %g secs)\n', NumChosen, cputime-T);
   end
end



