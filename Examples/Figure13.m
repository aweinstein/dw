% Figure13.m
%
% Generate Figure 13 in section 7.3.
%

close all;


if ~exist('SphereDataLoaded')
   fprintf('Figure13.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end

Coeffs = DWCoeffs(Tree, G);
Best   = DWBest(Coeffs);
BestList = DWUnpack(Best);

Wavelet = DWWavelet(Coeffs);
WaveletList = DWUnpack(Wavelet);
EigCoeffs = Eigs*G;

BestCoeffs = flipud(sort(abs(BestList(:,4))));
WaveletCoeffs = flipud(sort(abs(WaveletList(:,4))));
EigCoeffs = flipud(sort(abs(EigCoeffs)));


figure;
plot(log10(abs(BestCoeffs)));
hold on;
plot(log10(abs(WaveletCoeffs)), 'r');

figure;
plot(BestCoeffs(1:50));
hold on;
plot(WaveletCoeffs(1:50), 'r');

%figure;
%plot(BestCoeffs(1:200));
%hold on;
%plot(EigCoeffs(1:200), 'r');



fprintf('Figure 1: comparison of the top 50 best basis coefficients to the\n');
fprintf('top 50 wavelet basis coefficients.\n');

fprintf('Figure 2: log plot of the best basis coefficients versus the wavelet\n');
fprintf('basis coefficients.\n');


