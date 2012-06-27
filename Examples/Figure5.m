% Figure5.m
%
% Generate Figure 5 in section 5.1.
%


if ~exist('AnisotropicDataLoaded')
   fprintf('Figure4.m: Anistropic data not found; run "AnistropicData.m" first\n');
   return;
end

close all;

figure;
plot(Tree{7,2}.ExtBasis(:,1));
fprintf('figure 1: Anistropic wavelet function at position (7,2,1)\n');

figure;
plot(Tree{7,2}.ExtBasis(:,2));
fprintf('figure 2: Anistropic wavelet function at position (7,2,2)\n');

figure;
plot(Tree{7,3}.ExtBasis(:,1));
fprintf('figure 3: Anistropic wavelet packet function at position (7,3,1)\n');

figure;
plot(Tree{7,4}.ExtBasis(:,7));
fprintf('figure 4: Anistropic wavelet packet function at position (7,4,7)\n');


% print the caption
fprintf('\nCaption: \n\n');
fprintf('Some level 7 wavelet and wavelet packet functions on the anistropic circle.  The\n');
fprintf('indices are, from left to right, (7,2,1), (7,2,2), (7,3,1) and (7,4,7).\n');



