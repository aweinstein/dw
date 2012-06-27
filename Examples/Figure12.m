% Figure12.m
%
% Generate Figure 12 in section 7.3.
%

if ~exist('SphereDataLoaded')
   fprintf('Figure12.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end


figure;
h = axes;
DrawSphereFcn(Points, G);
set(h, 'View', [30 60])

figure;
h = axes;
DrawSphereFcn(Points, G);
set(h, 'View', [171.5 28])

fprintf('\nCaption:\n\n');
fprintf('Two different views of the function G on the sphere.\n');
