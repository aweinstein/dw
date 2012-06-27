% Figure14.m
%
% Generate Figure 14 in section 7.3.
%

if ~exist('SphereDataLoaded')
   fprintf('Figure14.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end

M = 200;

Coeffs = DWCoeffs(Tree, G);
Best   = DWBest(Coeffs);
BestList = DWUnpack(Best);

BestM = BestList(1:M,:);
Best = DWPack(Tree, BestM);

ReconG1 = DWRecon(Tree, Best);
figure;
DrawSphereFcn(Points, ReconG1);

EigCoeffs = Eigs*G;
[junk idxs] = sort(abs(EigCoeffs),'descend');
EigCoeffs(idxs(M+1:N)) = 0.0;
ReconG2 = Eigs'*EigCoeffs;


figure;
DrawSphereFcn(Points, ReconG2);

fprintf('Figure1: Reconstruction of G from the top 200 best basis coefficients.\n');
fprintf('Figure2: Reconstruction of G from the top 200 eigenfunction coefficients.\n');

fprintf('\n');
fprintf('Difference in norm between the best basis reconstructed G and G = %g\n', norm(G-ReconG1));
fprintf('Difference in norm between the eigenfunction reconstructed G and G = %g\n', norm(G-ReconG2));
