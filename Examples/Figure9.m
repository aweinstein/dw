% Figure9.m
%
% Generate Figure 9 in section 7.2.
%

if ~exist('SphereDataLoaded')
   fprintf('Figure6.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end


figure;
h = axes;
DrawSphereFcn(Points, F);
set(h, 'View', [30 60])

figure;
h = axes;
DrawSphereFcn(Points, F);
set(h, 'View', [171.5 28])

%fprintf('\nCaption:\n\n');
%fprintf('Two different views of the function F on the sphere.\n);