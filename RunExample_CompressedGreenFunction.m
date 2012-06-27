% Create diffusion operator
fprintf('Creating diffusion operator T...\n');
T=MakeCircleDiffusion(256);              % Create the diffusion operator
T(5:10,10:20) = 0.1;                    % Make the diffusion not symmetric

% Compute the top eigenvector, which is going to be the kernel of 1-T
fprintf('Computing kernel of 1-T...\n');
[lV lS]=eigs(T,1);
T = T / lS;

% Create the diffusion wavelets
fprintf('Creating diffusion wavelets...\n');
Tree = DWPTree(T,20,eps,struct('Wavelets',false,'ExtendBases',false,'StopFcns',1));

% Pick a vector F for which we'll compute (I-T)^{-1}F
F=[zeros(1,floor(size(T,1)/2)),floor(size(T,1)/2):-1:1]'; figure;plot(F);title('Vector to invert');

% Compute the inverse using the compressed Green's function
F_inv = DWApplyCompressedGreenFunction( Tree, T, F ,[], [], 15 );    % Compute the inverse

% Plot results
figure;plot(F_inv);title('Inverse vector');
figure;plot(F);hold on;plot((speye(size(T))-T)*(F_inv),'r');title('Original vector and image of the inverse');
figure;plot(log10(abs(F-((speye(size(T))-T)*(F_inv)))));title('Difference between original and image of the inverse');
lDiffV=F-((speye(size(T))-T)*(F_inv));
fprintf('Norm of error outside the null space:%f',norm(lDiffV-(lDiffV'*lV)*lV)),               % This is the norm after projecting out the null space
