function vX = DWApplyCompressedGreenFunction( cTree, cT, cF, cGamma, cD, cMaxPower )

%
% function vX = DWApplyCompressedGreenFunction( cTree, cT, cF, cGamma, cD, cMaxPower )
%
% IN:
%   cTree   : wavelet tree
%   cF      : function as N column vector to which to apply the Green's function 
%               (I-cGamma*diag(cD)^{-1/2}*cT*diag(cD)^{1/2})^{-1}, where T was used to build cTree
%             cF must be in the range of (I-cGamma*diag(cD)^{-1/2}*T*diag(cD)^{1/2}). If not, a solution does not exist, and the vector vX returned
%             satisfies
%                (I-cGamma*diag(cD)^{-1/2}*cT*diag(cD)^{1/2})vX-vF \in kernel((I-cGamma*diag(cD)^{-1/2}*T*diag(cD)^{1/2}))
%   [cGamma]: constant. Default 1.
%   [cD]    : N column vector.
%   [cMaxPower] : 2^cMaxPower is the highest power used, as in Schultz' formula. Default: size(cTree,1). TODO: this should be adaptive and sensible.
%
% OUT:
%   vX      : column vector of result of applying Green's function to cF, i.e. vX solves
%                   (I-cGamma*diag(cD)^{-1/2}*cT*diag(cD)^{1/2})vX = vF
%             If cF is not in the range of Green's operator, the difference between the left- and right-hand sides is in the kernel.
%
%
% EXAMPLE:
%   cDiffusionMatrix=MakeCircleDiffusion(64);              % Create the diffusion operator
%   cDiffusionMatrix(5:10,10:20) = 0.1;                    % Make the diffusion not symmetric
%   [lV lS]=eigs(cDiffusionMatrix,1);
%   cDiffusionMatrix = cDiffusionMatrix / lS;
%   vTree = DWPTree(cDiffusionMatrix,20,eps,struct('Wavelets',false,'ExtendBases',false,'StopFcns',1));
%   lVector=[zeros(1,floor(size(cDiffusionMatrix,1)/2)),floor(size(cDiffusionMatrix,1)/2):-1:1]'; figure;plot(lVector);title('Vector to invert');
%   vX = DWApplyCompressedGreenFunction( vTree, cDiffusionMatrix, lVector ,[], [], 15 );    % Compute the inverse
%   figure;plot(vX);title('Inverse vector');
%   figure;plot(lVector);hold on;plot((speye(size(cDiffusionMatrix))-cDiffusionMatrix)*(vX),'r');title('Original vector and image of the inverse');
%   figure;plot(lVector-((speye(size(cDiffusionMatrix))-cDiffusionMatrix)*(vX)));title('Difference between original and image of the inverse');
%   lDiffV=lVector-((speye(size(cDiffusionMatrix))-cDiffusionMatrix)*(vX));
%   norm(lDiffV-(lDiffV'*lV)*lV),               % This is the norm after projecting out the null space
%   
%
% SC:
%   MM  :   5/23/2006   : New version for new DW code
%   MM  :   4/11/2007   : New version for new(er) DW code
%
% TODO:
%   One can work all the time in wavelet domain and do only one inversion at the end: need to modify DWApplyDytadicPower and the loop below.
%   This will improve the constants of the algorithm quite dramatically.
%
%
% (c) Copyright Yale University, 2006, Mauro Maggioni
% (c) Copyright Duke University, 2007, Mauro Maggioni
%

lN = length(cF);

if (nargin<3) | (isempty(cGamma)),
    cGamma = 1;
end;
if (nargin<4) | (isempty(cD)),
    cD = ones(lN,1);
end;
if (nargin<5) | (isempty(cMaxPower)),
    cMaxPower = size(cTree,1);
end;

% Initialize result
vX = cF;

lSqrtD      = sqrt(cD);
lInvSqrtD   = lSqrtD.^(-1);
lSqrtD      = spdiags(lSqrtD,0,lN,lN);
lInvSqrtD   = spdiags(lInvSqrtD,0,lN,lN);

% Go through all levels
for lk = 0:1:cMaxPower,
    lTemp = DWApplyDyadicPower(cTree,lSqrtD*vX,lk, cT);
    if length(lTemp)==0,
        break;
    end;
    vX = vX+cGamma^(2^lk)*lInvSqrtD*lTemp;    
end;

% Subtract from cF the component in the kernel of the operator
%[lV lS]=eigs(cT,1,'LM',struct('disp',0));
%vX = vX-(vX'*lV)*lV;

return;