% Figure3.m
%
% Generate Figure 3 in section 5.1.
%

if ~exist('AnisotropicDataLoaded')
   fprintf('Figure4.m: Anistropic data not found; run "AnistropicData.m" first\n');
   return;
end

close all;
figure;

plot(Impedance);

fprintf('Caption\n\n');
fprintf('Impedance of the anistropic diffusion operator T: high on one part of the\n');
fprintf('circle and low on another.\n');


