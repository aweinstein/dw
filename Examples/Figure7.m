% Figure7.m
%
% Generate Figure 7 in section 7.1.
%

if ~exist('AnisotropicDataLoaded')
   fprintf('Figure4.m: Anistropic data not found; run "AnistropicData.m" first\n');
   return;
end

close all;

figure;
plot(F);
fprintf('figure 1: the function F\n');

figure;
plot(flipud(F));
fprintf('figure 2: the reflection of F\n');

fprintf('\nCaption:\n\n');
fprintf('The function F (on the left) and its reflection.\n');





