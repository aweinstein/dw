% AnisotropicData.m
%
% Script for generating the anistropic diffusion data.  This data is used in
% sections 5.1 and 7.1.  Once the data has been generated it is saved to the
% file 'anistropicxxx.mat' where xxx is the number of points used.
%
%
% Variables:
%    N         = number of points on the circle
%    T         = anisotropic diffusion
%    Tree      = diffusion wavelet packet tree
%    Impedance = impedance of the operator
%    F         = function for compression example


clear;      % get rid of troublesome data
N = 512;    % number of points to use

filename = sprintf('anisotropic%d.mat', N);

%{
% if the file anistropic.mat exists, just load the data from there
%}
if exist(filename)
   fprintf('Anisotropic.m: loading anisotropic data ... '); time=cputime;
   load(filename);
   fprintf('%g seconds\n', cputime-time);
else
   fprintf('Anisotropic.m: could not find %s, generating data from scratch ... \n', filename);
   %{
   % generate vector of impedances
   %}
   lNumberOfPoints = N;
   lImpedance =0.5+0.49*sin(2*pi/(lNumberOfPoints-1)*(1:lNumberOfPoints));
   lImpedance = lImpedance * 10;

   % allocate memory for the diffusion matrix
   cDiffusionMatrix = zeros(lNumberOfPoints,lNumberOfPoints);

   % create the diffusion matrix by hand
   cDiffusionMatrix(1,[1 2 lNumberOfPoints]) = [4-(lImpedance(1)+lImpedance(lNumberOfPoints)) lImpedance(1) lImpedance(lNumberOfPoints)];      %-(lImpedance(1)+lImpedance(lNumberOfPoints))
   for lk = 2:lNumberOfPoints-1
      cDiffusionMatrix(lk,lk-1:lk+1) = [lImpedance(lk-1) 4-(lImpedance(lk-1)+lImpedance(lk)) lImpedance(lk)];       %-(lImpedance(lk-1)+lImpedance(lk))
   end;
   cDiffusionMatrix(lNumberOfPoints,[1 lNumberOfPoints-1 lNumberOfPoints]) = [lImpedance(lNumberOfPoints) lImpedance(lNumberOfPoints-1) 4-(lImpedance(lNumberOfPoints)+lImpedance(lNumberOfPoints-1))]; %-(lImpedance(lNumberOfPoints)+lImpedance(lNumberOfPoints-1))

   % normalize and generate the set of points
   T = Bimarkov(sparse(cDiffusionMatrix));%sparse(GraphNormalize(cDiffusionMatrix,1e-10));
   Thetas = 0:2*pi/N:2*pi-2*pi/N;
   Points    = [cos(Thetas) sin(Thetas)];

   % construct the tree
   Tree = DWPTree(T, 8, 1e-5, struct('StopFcns', 10));

   % construct the function used for the compression example
   center1 = N*1/4;
   center2 = N*3/4;
   sigma   = 850;

   Bump1 = exp(-1/sigma*((1:N)-center1).^2);
   Bump2 = exp(-1/sigma*((1:N)-center2).^2);

   % generate the function F used for the compression example
   F = (sin(Thetas.^2).*Bump2+sin(17*Thetas.^5).*Bump1)';


   % indicator that we have the anisotropic data loaded in memory
   AnisotropicDataLoaded = 1;
   Impedance = lImpedance;

   % save the data to a file
   save(filename);
end
