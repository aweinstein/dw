function vF = DWApplyDyadicPower( cTree, cF, cPower, cT )

%
% function vF = DWApplyDyadicPower( cTree, cF, cPower )
%
% Computes a power
%
% IN:
%   cTree   :   Wavelet tree; should come from DWPTree, with expansion type at least 'Partial'
%   cF      :   function, in a column vector
%   cPower  :   will apply 2^cPower
%
% OUT:
%   vCoeffs :   scaling and wavelet transform coefficients of (T^{2^(cPower)})cF
%
% EXAMPLE:
%   cDiffusionMatrix=MakeCircleDiffusion(256);              % Create the diffusion operator
%   vTree = DWPTree(cDiffusionMatrix,20,eps,struct('Wavelets',false,'ExtendBases',1,'StopFcns',1));
%   lVector=[zeros(1,floor(size(cDiffusionMatrix,1)/2)),floor(size(cDiffusionMatrix,1)/2):-1:1]'; % Create the vector to take powers of
%   figure;plot(lVector);title('Vector to take powers of');
%   lPower = 15;
%   vF1 = DWApplyDyadicPower( vTree, lVector, lPower, cDiffusionMatrix );    % Compute a power
%   figure;plot(cDiffusionMatrix^(2^(lPower))*lVector,'b');hold on;plot(vF1,'r');
%   figure;plot(log10(abs(cDiffusionMatrix^(2^(lPower))*lVector-vF1)));
%
% SC:
%   MM  :   5/23/06     : update for the new DW code.
%   MM  :   4/8/07      : update for the new(er) DW code.
%
%
% (c) Copyright Yale University, Mauro Maggioni
% (c) Copyright Duke University, Mauro Maggioni
%

% All we need to do is to apply the product of Op_j (=M_j in the Diffusion Wavelet paper), j=0,\dots,cPower
lF = cF;

if cPower==0,
    vF = cT*lF;
elseif cPower==1,
    vF = cT*(cT*lF);
else
    lF = cT*lF;

    lLastLevel = min([size(cTree,1),cPower]);

    % Use the compressed form as much as possible
    for lk = 1:lLastLevel,
        lF = cTree{lk,1}.Op*lF;
        % Debug stuff
        %lk,norm(cTree{lk,1}.ExtBasis*lF-cT^(2^lk)*cF),figure(87);clf;plot(cTree{lk,1}.ExtBasis*lF);hold on;plot(cT^(2^lk)*cF,'r');pause;fprintf('Iteration %d, length %d \n',lk,length(lF));
    end;

    % If still need to take even greater power, use the highest power in compressed
    % form and iterate it (i.e. work in the top scaling function space.
    for lk = 1:(2^(cPower-size(cTree,1)+1)-2)
        lF = cTree{size(cTree,1),1}.T{1}*lF;
    end;

    % Return the results in form of scaling function coefficients, in tree format
    vCoeffs = zeros(size(lF,1),4);
    vCoeffs(:,1) = lLastLevel;              % Level
    vCoeffs(:,2) = 1;                       % Node (1=scaling functions)
    vCoeffs(:,3) = (1:size(vCoeffs,1))';    % All scaling function coefficients will be set
    vCoeffs(:,4) = lF;                      % These are the scaling function coefficients

    vF = DWRecon(cTree,DWPack(cTree,vCoeffs));        % Pack the coefficients in the tree structure and reconstruct
end;

return;