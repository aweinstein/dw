% Figure10.m
%
% Generate Figure 10 in section 7.2.
%

close all;

if ~exist('SphereDataLoaded')
   fprintf('Figure10.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end

Coeffs = DWCoeffs(Tree, F);
Best   = DWBest(Coeffs);
BestList = DWUnpack(Best);

Wavelet = DWWavelet(Coeffs);
WaveletList = DWUnpack(Wavelet);
EigCoeffs = Eigs*F;

BestCoeffs = flipud(sort(abs(BestList(:,4))));
WaveletCoeffs = flipud(sort(abs(WaveletList(:,4))));
EigCoeffs = flipud(sort(abs(EigCoeffs)));

figure;
plot(log10(abs(BestCoeffs)));
hold on;
plot(log10(abs(EigCoeffs)), 'r');

figure;
plot(log10(abs(BestCoeffs)));
hold on;
plot(log10(abs(WaveletCoeffs)), 'r');



fprintf('Figure 1: log plot of the decay of the magnitude of coefficients\n');
fprintf('in the best basis (blue) and in the basis of eigenfunctions(red).\n');

fprintf('Figure 2: log plot of the decay of the magnitude of coefficients\n');
fprintf('in the best basis (blue) and in the basis of wavelets (red).\n');

