% SphereData.m
%
% Script to generate the data for sphere examples discussed in sections
% 5.2, 7.2, 7.3, and 7.4.
%
% Variables:
%
%    Points = set of 2000 points on the sphere
%    T      = diffusion operator on Points
%    Eigs   = basis of oversampled eigenfunctions
%    Tree   = 14 level diffusion wavelet packet tree for the operator T
%    F      = first compression example function
%    G      = second compression example function

clear;      % we'll have more than enough data as it is
N = 2000;   % number of points to use
M = 10000;  % number of points to use for oversampling

filename = sprintf('sphere%d.mat', N);

% check to see if the data has already been generated
if exist(filename)
    fprintf('SphereData.m: loading %s data file ...', filename);
    load(filename);
    fprintf('\n');
else
    
    fprintf('SphereData.m: could not find %s, generating data from scratch ... \n', filename);
    
    Points = GeneratePoints('Sphere', N);
    OverPoints = [Points, GeneratePoints('Sphere', M-N)];   % oversampled version
    
    % construct a diffusion operator and an oversampled diffusion operator
    T = MakeDiffusion(Points, 10, struct('Delta', 160, 'Threshold', 1e-3, 'Normalization', 'bimarkov'));
    OverT = MakeDiffusion(OverPoints, 10, struct('Delta', 160, 'Threshold', 1e-3, 'Normalization', 'bimarkov'));
    
    %    T = T^2;
    %    OverT = OverT^2;
    
    % compute eigenvalues and vectors for the oversampled sphere
    fprintf('SphereData.m: solving eigenproblem for overampled T ... '); TIME=cputime;    
    tic;    
    [V, D] = eigs(OverT,N);
    % flip the eigenvectors
    V = fliplr(V);        
    fprintf('%g seconds\n', cputime-TIME);
    
    % project the eigenfunctions down
    ProjectedV = zeros(N, N); warning off;
    for j=1:N
        % vectorize call to griddata ?????
        ProjectedV(:, j) = griddata3(OverPoints(1,:), OverPoints(2,:), OverPoints(3,:), V(:,j), Points(1,:), Points(2,:), Points(3,:), 'nearest');
    end
    
    % now re-orthogonalize the projected vectors
    fprintf('SphereData.m: reorthogonalizing projected eigenvectors ... ');
    %Eigs = gs(ProjectedV);
    SpQR = gsqr(ProjectedV, struct('StopN', N, 'StopPrecision', -1.0));
    Idxs = SpQR.Idxs;
    DIdxs = SpQR.DIdxs;
    Eigs = V;
    
    % now compute the tree
    Tree = DWPTree(T, 12, 1e-5, struct('StopFcns', 4, 'ExtendBases', false));
    
    % build the two functions used in the compression examples
    Points=Points';
    [Theta Phi Rho] = cart2sph(Points(:,1), Points(:,2), Points(:,3));
    
    % F is trigonemtric poly + exponential at different frequncies
    %   F = sin(3*Theta.^8) .* exp(-5*Phi.^2);
    % build F
    F = sin(Theta);	% background
    
    % add swirl
    Pt = Points(26,:);
    
    atria = nn_prepare(Points);
    [count neighbors] = range_search(Points, atria, Pt, 1.0);
    
    Delta = 10;
    idxs = neighbors{1};
    dists = neighbors{2};
    
    ips = Points(idxs, :)*Pt';
    
    %F1(neighbors{1}(j)) = Points(neighbors{1}(j),1)^4*sin(20*ip) * exp(-Delta*neighbors{2}(j)^2);
    F(idxs) = F(idxs)+Points(idxs,1).^4.*sin(20*ips).*exp(-Delta*dists.^2)';
    
    idxs = find(0 < Theta & Theta < pi/4 & 0 < Phi & Phi < pi/6);
    F(idxs) = -1.0;
    
    G = zeros(N,1);
    G = G+.5*cos(.6*Theta).*exp(-3*Phi);
    
    %G = cos(Phi).^3-sin(Theta.^2)+cos(3*Phi-Theta);
    %idxs = find(-pi/2 <= Phi & Phi <= pi/2 & pi/4<= Theta & Theta <= pi/2);
    %G(idxs) = 1.4;
    
    fprintf('\n\nSphereData.m: saving sphere data to %s ...', filename); TIME=cputime;
    SphereDataLoaded = 1;
    save(filename);
    fprintf('%g seconds\n', cputime-TIME);
end

