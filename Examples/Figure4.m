% Figure4.m
%
% Generate Figure 4 in section 5.1.
%


if ~exist('AnisotropicDataLoaded')
   fprintf('Figure4.m: Anistropic data not found; run "AnisotropicData.m" first\n');
   return;
end

close all;

figure;
plot(Tree{7,1}.ExtBasis(:,1));
fprintf('figure 1: Anistropic scaling function at position (7,1,1)\n');

figure;
plot(Tree{7,1}.ExtBasis(:,5));
fprintf('figure 2: Anistropic scaling function at position (7,1,5)\n');

figure;
plot(Tree{7,1}.ExtBasis(:,111));
fprintf('figure 3: Anistropic scaling function at position (7,1,111)\n');

figure;
plot(Tree{7,1}.ExtBasis(:,128));
fprintf('figure 4: Anistropic scaling function at position (7,1,1)\n');


% print the caption
fprintf('\nCaption: \n\n');
fprintf('Some level 7 scaling functions on the anistropic circle.  Notice that the\n');
fprintf('functions exhibit high frequency behavior in the high impedance region and low\n');
fprintf('frequency behavior in the low impedance region.  The functions are, from left to right\n');
fprintf('at indices (7,1,1), (7,1,5), (7,1,111) and (7,1,1).\n');
