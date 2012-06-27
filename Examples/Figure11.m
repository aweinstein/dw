% Figure11.m
%
% Generate Figure 11 in section 7.2.
%

if ~exist('SphereDataLoaded')
   fprintf('Figure6.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end

M = 200;


Coeffs = DWCoeffs(Tree, F);
Best   = DWBest(Coeffs);
BestList = DWUnpack(Best);

BestM = BestList(1:M,:);
Best = DWPack(Tree, BestM);

ReconF1 = DWRecon(Tree, Best);
figure;
DrawSphereFcn(Points, ReconF1);

EigCoeffs = Eigs*F;
[junk idxs] = sort(abs(EigCoeffs),'descend');
EigCoeffs(idxs(M+1:N)) = 0.0;
ReconF2 = Eigs'*EigCoeffs;


figure;
DrawSphereFcn(Points, ReconF2);

fprintf('Figure1: Reconstruction of F from the top 200 best basis coefficients.\n');
fprintf('Figure2: Reconstruction of F from the top 200 eigenfunction coefficients.\n');

fprintf('\n');
fprintf('Difference in norm between the best basis reconstructed F and F = %g\n', norm(F-ReconF1));
fprintf('Difference in norm between the eigenfunction reconstructed F and F = %g\n', norm(F-ReconF2));
