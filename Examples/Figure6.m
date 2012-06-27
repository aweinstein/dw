% Figure6.m
%
% Generate Figure 6 in section 5.2.
%

if ~exist('SphereDataLoaded')
   fprintf('Figure6.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end

close all;
figure;
DrawSphereFcn(Points, DWBasisFcn(Tree, 3, 2, 3));
figure;
DrawSphereFcn(Points, DWBasisFcn(Tree, 9, 2, 5));
figure;
DrawSphereFcn(Points, DWBasisFcn(Tree, 5, 3, 8));
figure;
DrawSphereFcn(Points, DWBasisFcn(Tree, 9, 5, 6));

fprintf('Figure 1: diffusion wavelet with index (3,2,3)\n');
fprintf('Figure 2: diffusion wavelet with index (9,2,5)\n');
fprintf('Figure 3: diffusion wavelet packet with index (5,3,8)\n');
fprintf('Figure 4: diffusion wavelet with index (9,5,6)\n');



fprintf('\nCaption:\n\n');
fprintf('Diffusion wavelets and wavelet packets on the sphere.  The index of\n');
fprintf('each function are: upper left (3,2,3), upper right (9,2,5), \n');
fprintf('lower left (5,3,8), lower right (9,5,6)\n');
