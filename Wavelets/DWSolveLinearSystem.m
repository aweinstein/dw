function vX = DWSolveLinearSystem( cA, cb, cOpts )

%
%

% EXAMPLE:
%   cDiffusionMatrix=CreateCircleDiffusion(64);              % Create the diffusion operator
%   cDiffusionMatrix(5:10,10:20) = 0.1;
%   [lV lS]=eigs(cDiffusionMatrix,1);
%   cDiffusionMatrix = cDiffusionMatrix / lS;
%   lVector=[zeros(1,floor(size(cDiffusionMatrix,1)/2)),floor(size(cDiffusionMatrix,1)/2):-1:1]'; figure;plot(lVector);title('Vector to invert');
%   vX = DWSolveLinearSystem( cDiffusionMatrix, lVector, struct('Type','Normal') );
%   figure;plot(vX);title('Inverse vector');
%   A = speye(size(cDiffusionMatrix))-cDiffusionMatrix;
%   figure;plot(lVector);hold on;plot(A*vX,'r');title('Original vector and image of the inverse');
%   figure;plot(lVector-A*vX);title('Difference between original and image of the inverse');
%   lSol=pinv(full(A))*lVector;
%   norm(lSol-vX),norm(lSol),norm(vX),norm(A*lSol-lVector)/norm(lVector),norm(A*vX-lVector)/norm(lVector)           
%

if nargin<3,
    cOpts = [];
end;

if ~isfield(cOpts,'Type'),
    cOpts.Type = 'None';
end;
if ~isfield(cOpts,'Levels'),
    cOpts.Levels = 20;
end;
if ~isfield(cOpts,'Powers'),
    cOpts.Powers = [];
end;
if ~isfield(cOpts,'Precision'),
    cOpts.Precision = eps;
end;
if ~isfield(cOpts,'Threshold'),
    cOpts.Threshold = 0;
end;
if ~isfield(cOpts,'IPThreshold'),
    cOpts.IPThreshold = 0;
end;
if ~isfield(cOpts,'Reorthogonalize'),
    cOpts.Reorthogonalize = 1;
end;

if strcmpi(cOpts.Type,'None'),
    lT = cA;
    lF = cb;
elseif strcmpi(cOpts.Type,'Normal'),
    % Build the normal equations
    lAt = cA';
    lT = cA+lAt-lAt*cA;
    lF = cb-lAt*cb;
    clear lAt;
end;

% Make sure lT has operator norm =1
lS = svd(full(lT));
lT = lT/lS(1);

% Create the diffusion wavelet tree
vTree = DWPTree(lT,cOpts.Levels,cOpts.Precision,struct('Wavelets',false,'ExtendBases',false,'StopFcns',1,'GSOptions',struct('Threshold',cOpts.Threshold,'IPThreshold',cOpts.IPThreshold,'Reorthogonalize',cOpts.Reorthogonalize,'Quiet',1)));

% Invert the linear system
vX = DWApplyCompressedGreenFunction( vTree, lT, lF, lS(1), [], cOpts.Powers );

return;