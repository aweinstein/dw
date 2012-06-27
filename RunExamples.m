%
% This file runs several examples of diffusion maps
%
clear all;close all;clc

global X Data;

pExampleNames   = { 'Circle', ...
    '2-d sphere', ...
    'Sphere and segment', ...
    'Spiral and plane', ...
    'Meyerstaircase', ...
    'D-Cube with low sampling rate', ...
    'S-Manifold', ...
    '10-d cube with high noise', ...
    '9-d sphere with high noise', ...
    'Two lines and a plane', ...
    'Isomap faces', ...
    'CBCL faces, I', ...
    'CBCL faces, II', ...
    'Science News articles', ...
    'MNIST digits', ...
    'Patches from Lena Image', ...
    'Anisotropic Diffusion on Circle'};

fprintf('\n\n Select example to run:\n');
for k = 1:length(pExampleNames),
    fprintf('\n [%d] %s',k,pExampleNames{k});
end;
fprintf('\n\n  ');

while true,
    if (~exist('pExampleIdx') || isempty(pExampleIdx) || pExampleIdx==0),
        try
            pExampleIdx = input('');
            pExampleIdx = str2num(pExampleIdx);
        catch
        end;
    end;
    if (pExampleIdx>=1) && (pExampleIdx<=length(pExampleNames)),
        break;
    else
        fprintf('\n %d is not a valid Example. Please select a valid Example above.',pExampleIdx);
        pExampleIdx=0;
    end;
end;

%% Set parameters for constructing the graph
EstDimOpts = struct('NumberOfTrials',15,'verbose',0,'MAXDIM',100,'MAXAMBDIM',100,'Ptwise',false,'NetsOpts',[],'UseSmoothedS',false, 'EnlargeScales',true );

%% Set parameters for data generation and generates the data set
switch pExampleIdx
    case 1
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',1,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1/sqrt(100));
    case 2
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',2,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
    case 3
        XName = 'SphereAndLine'; XNickName = 'Sphere and Line';
        XOpts = struct('Dim','2.5','NoiseType','Gaussian','NoiseParam',0.00);
        EstDimOpts.Ptwise = true;
    case 4
        XName = 'SpiralAndPlane'; XNickName = 'Spiral and Plane';
        XOpts = struct('NumberOfPoints',1100,'Dim','1.5');
        EstDimOpts.Ptwise = true;
        Labels=[ones(300,1);2*ones(800,1)];
    case 5
        XName = 'MeyerStaircase'; XNickName = 'Z';
        XOpts = struct('NumberOfPoints',500,'Dim',1000,'MeyerStepWidth',20,'EmbedDim',1000,'NoiseType','Gaussian','NoiseParam',0.05/sqrt(1000));
    case 6
        XName = 'D-Cube'; XNickName = 'Q';
        XOpts = struct('NumberOfPoints',250,'Dim',6,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01/sqrt(100));
    case 7
        XName = 'S-Manifold'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',2,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01);
    case 8
        XName = 'D-Cube'; XNickName = 'Q';
        XOpts = struct('NumberOfPoints',1000,'Dim',10,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
    case 9
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',1000,'Dim',9,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
        EstDimOpts.EnlargeScales = true;
    case 10
        XName = 'TwoLinesAndAPlane'; XNickName = 'Two Lines and a Plane';
        XOpts = struct('NumberOfPoints',800,'Dim','1.5');
        EstDimOpts.Ptwise = true;
        Labels= [ones(400,1);2*ones(400,1)];
    case 11
        XName = 'Isomapfaces'; XNickName = 'Isomap Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 12
        XName = 'CBCLFaces1'; XNickName = 'CBCLFaces1 Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 13
        XName = 'CBCLFaces2'; XNickName = 'CBCLFaces2 Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 14
        XName = 'ScienceNews'; XNickName = 'ScienceNews Articles';
        XOpts = struct(); %struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
    case 15
        XName = 'BMark_MNIST'; XNickName = 'MNIST Digits';
        XOpts = struct('NumberOfPoints',200,'MnistOpts',struct('Sampling', 'RandN', 'QueryDigits',0:9, 'ReturnForm', 'vector'));
    case 16
        XName = 'ImagePatches'; XNickName = 'Lena Patches';
        XOpts = struct('ImageFileName','Lena.jpg','PatchSize',16,'DownsampleImage',8);
    case 17
        [ Points, Basis, Diffusion, Tree ] = MakeExample( 'Anisotropic',1 );
        return;
end

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
    title(sprintf('Diffusion scaling function at scale %d',k));
end;

figure;
for j = 1:4,
    subplot(2,2,j);
    scatter3(G.EigenVecs(:,2),G.EigenVecs(:,3),G.EigenVecs(:,4),15,Tree{j,2}.ExtBasis(:,1));hold on;
    title(sprintf('Diffusion wavelet at scale %d',k));
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