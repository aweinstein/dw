% Figure 8
%
% Generate figure 8 in section 7.1


if ~exist('AnisotropicDataLoaded')
   fprintf('Figure4.m: Anistropic data not found; run "AnistropicData.m" first\n');
   return;
end

close all;

% find its best basis
Coeffs = DWCoeffs(Tree, F);
Basis  = DWBest(Coeffs);
list   = DWUnpack(Basis);


figure;
plot(abs(list(1:50,4)));
hold on;

Coeffs = DWCoeffs(Tree, flipud(F));
Basis  = DWBest(Coeffs);
list   = DWUnpack(Basis);

plot(abs(list(1:50,4)), 'r');

fprintf('\nCaption:\n\n');
fprintf('Comparison of the magnitude of the largest 50 coefficients of F and\n');
fprintf('its reflection in their best wavelet packet bases.\n');


