EstDimOpts = struct('NumberOfTrials', 15, 'verbose', 0, 'MAXDIM', 100, ...
                    'MAXAMBDIM', 100, 'Ptwise', false, 'NetsOpts', [], ...
                    'UseSmoothedS', false, 'EnlargeScales', true);
XName = 'D-Sphere'; XNickName = 'S';
XOpts = struct('NumberOfPoints', 1000, 'Dim', 1, 'EmbedDim', 100, ...
               'NoiseType', 'Gaussian', 'NoiseParam', 0.1/sqrt(100));

%% Generate the data set
fprintf('\nGenerating %s data...', XName);
[X,GraphDiffOpts] = GenerateDataSets( XName, XOpts);
X = 1*X;
fprintf('done.');

%% Compute Graph Diffusion, and display diffusion embedding
fprintf('\n\nConstructing graph and diffusion map...\n');
G = GraphDiffusion(X, 0, GraphDiffOpts);                  % These are for debugging purposes
fprintf('done.\n');

figure;
subplot(1,2,1);plot3(G.EigenVecs(:,2),G.EigenVecs(:,3),G.EigenVecs(:,4),'.');title('Diffusion embedding (2,3,4)');
subplot(1,2,2);plot3(G.EigenVecs(:,5),G.EigenVecs(:,6),G.EigenVecs(:,7),'.');title('Diffusion embedding (5,6,7)');

%% Compute Diffusion Wavelets
Tree = DWPTree (G.T^4, 12, 1e-4, struct('Wavelets',true,'OpThreshold',1e-2,'GSOptions',struct('StopDensity',1,'Threshold',1e-2)));

%% Show some diffusion scaling functions, in regular scale and log10 scale to show exp. decay
figure;
for j = 1:4,
    subplot(2,2,j);
    scatter3(G.EigenVecs(:,2),G.EigenVecs(:,3),G.EigenVecs(:,4),15,Tree{j,1}.ExtBasis(:,1));hold on;
    title(sprintf('Diffusion scaling function at scale %d',j));
end;

figure;
for j = 1:4,
    subplot(2,2,j);
    scatter3(G.EigenVecs(:,2),G.EigenVecs(:,3),G.EigenVecs(:,4),15,Tree{j,2}.ExtBasis(:,1));hold on;
    title(sprintf('Diffusion wavelet at scale %d',j));
end;

% Show some compressed powers of T
figure;
for j = 1:4,
    subplot(2,4,j);
    imagesc(Tree{j,1}.T{1});title(sprintf('[T]_%d',j));colorbar;
    subplot(2,4,j+4);
    imagesc(log10(abs(Tree{j,1}.T{1})));title(sprintf('log_{10}(|[T]_%d|)',j));colorbar;
end;

fprintf('\n');

return;

fprintf('\n');