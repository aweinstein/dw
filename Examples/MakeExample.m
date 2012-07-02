function [ Points, Basis, Diffusion, Tree ] = MakeExample( ExampleName, Force)

% function [ Points, Basis, Diffusion, Tree ] = MakeExample( ExampleName, Force )
%
% Create the data for one of the examples in the diffusion wavelet packet paper.
% This function first looks for a .mat file with the same name as the
% example and if found, loads the data from there.
%
% In:
%    ExampleName = String indicating which of the following examples to generate
%                  data for:
%
%                  'Circle'      = standard diffusion on the circle
%                  'Anisotropic' = an anistropic diffusion on the circle
%                  'Sphere'      =
%
%    Force       = If true, regenerate the data whether it is already
%                  stored or not.
%
% Uses:
%    DWPTree, GraphDiffusion
%
% Out:
%
%    Points      = An NxM matrix of N points on the graph.
%    Basis       = The delta basis.
%    Diffusion   = The diffusion operator on the delta basis.
%    Tree        = The diffusion wavelet tree.
%
%
% (c) Copyright Duke University
%

P = [];
B = [];
T = [];
Tree = [];

DebugPlots = 1;
EigPlot    = 1;

if nargin < 2
    Force = 1;
end


Options.Extend    = 1;
Options.DebugPlot = 0;
Options.DebugOut  = 1;
%Options.SplitFcn  = 'WPSplit';

if strcmpi(ExampleName, 'Circle')
    
    if ~exist('Circle.mat') || Force
        N          = 256;
        Delta      = 70;
        TPrecision = 1e-3;
        Levels     = 15;
        Precision  = 1e-10;
        
        fprintf('[-] Generating diffusion operator ...\n\n');
        
        % make some points
        Thetas = linspace(0,2*pi, N+1);
        Thetas = Thetas(1:N);
        
        Points    = [cos(Thetas);sin(Thetas)];
        Basis     = speye(N);
        
        Options.kNN=2;
        Options.Delta=70;
        Diffusion = GraphDiffusion(Points, Options.Delta, Options);
        
        %      Diffusion = MakeDiffusion(Points, Delta, TPrecision);
        
        if DebugPlots
            figure; scatter(Points(1,:), Points(2,:)); title('Circle: Points');
            figure; imagesc(Diffusion.T); title('Circle: image');
            
            if EigPlot
                %figure; plot(flipud(sort(abs(eigs(Diffusion.T, N-2))))); title('Circle: Eigs');
                figure; plot(flipud(sort(eigs(Diffusion.T, N-2)))); title('Circle: Eigs');
            end
            [V, D] = eigs(Diffusion.T, 4);
            figure; scatter(V(:,2), V(:,3)); title('Circle: Coifman-Lafon embedding');
        end
        
        Tree = [];
        
        fprintf('[-] Making diffusion wavelet tree ...\n\n');
        % generate the tree
        Tree = DWPTree(Diffusion.T,Levels,Precision,Options);
        
        % save the data to a mat file
        save 'Circle.mat' Points Basis Diffusion Tree N Precision Levels;
    else
        load 'Circle.mat';
    end
elseif strcmpi(ExampleName, 'Anisotropic')
    
    if ~exist('Anisotropic.mat') | Force
        
        N          = 256;
        Levels     = 10;
        
        %   TPrecision = 1e-3;
        Precision  = 1e-10;
        
        lNumberOfPoints = N;
        %lImpedance = 0.25+(1.9*((1:lNumberOfPoints)-128)/(lNumberOfPoints)).^6.*sin(pi/lNumberOfPoints*(1:lNumberOfPoints));
        lImpedance =0.5+0.49*sin(2*pi/(lNumberOfPoints-1)*(1:lNumberOfPoints));
        
        lImpedance = lImpedance * 10;
        %lImpedance(1:floor(lNumberOfPoints/2))=1;
        %lImpedance(floor(lNumberOfPoints/2)+1:lNumberOfPoints)=2;
        %lImpedance=0.5*ones(1,lNumberOfPoints);
        
        % Allocate memory for the diffusion matrix
        cDiffusionMatrix = zeros(lNumberOfPoints,lNumberOfPoints);
        
        cDiffusionMatrix(1,[1 2 lNumberOfPoints]) = [4-(lImpedance(1)+lImpedance(lNumberOfPoints)) lImpedance(1) lImpedance(lNumberOfPoints)];      %-(lImpedance(1)+lImpedance(lNumberOfPoints))
        for lk = 2:lNumberOfPoints-1
            cDiffusionMatrix(lk,lk-1:lk+1) = [lImpedance(lk-1) 4-(lImpedance(lk-1)+lImpedance(lk)) lImpedance(lk)];       %-(lImpedance(lk-1)+lImpedance(lk))
        end;
        cDiffusionMatrix(lNumberOfPoints,[1 lNumberOfPoints-1 lNumberOfPoints]) = [lImpedance(lNumberOfPoints) lImpedance(lNumberOfPoints-1) 4-(lImpedance(lNumberOfPoints)+lImpedance(lNumberOfPoints-1))]; %-(lImpedance(lNumberOfPoints)+lImpedance(lNumberOfPoints-1))
        
        % Normalize the diffusion matrix by making it row-stochastic
        %cDiffusionMatrix = cDiffusionMatrix * diag(sum(cDiffusionMatrix,2).^(-1),0);
        
        % Symmetrize the matrix
        %cDiffusionMatrix = cDiffusionMatrix*cDiffusionMatrix';
        
        Diffusion = sparse(GraphNormalize(cDiffusionMatrix,1e-10));
        Basis = speye(lNumberOfPoints);
        Thetas = linspace(0,2*pi, N+1);
        Thetas = Thetas(1:N);
        Points    = [cos(Thetas);sin(Thetas)];
        
        if DebugPlots
            figure; scatter(Points(1,:), Points(2,:)); title('Anisotropic Circle: Points');
            figure; imagesc(Diffusion); title('Anisotropic Circle: image');colorbar;
            if EigPlot
                figure; plot(flipud(sort(abs(eigs(Diffusion, N-2))))); title('Anisotropic Circle: Eigs');
            end
            
            [V, D] = eig(full(Diffusion));
            
            figure; scatter(V(:,end-1), V(:,end-2)); title('Anisotropic Circle: Coifman-Lafon embedding');
        end
        
        Tree = DWPTree(Diffusion,Levels,Precision,Options);
        
        save 'Anisotropic.mat' Points Basis Diffusion Tree N Precision Levels lImpedance;
    else
        load 'Anisotropic.mat';
    end    
else
    fprintf('I do not know how to make "%s.".\n', ExampleName);
    Points = [];
    Basis = [];
    Diffusion = [];
    Tree = {};
    return;
end


%% Show some diffusion scaling functions, in regular scale and log10 scale to show exp. decay
figure;
for j = 4:7,
    if ~isempty(Tree{j,1}.ExtBasis),
        subplot(2,2,j-3);
        plot(Points,Tree{j,1}.ExtBasis(:,5));hold on;
        title(sprintf('Diffusion scaling function at scale %d',j));
    end;
end;

figure;
for j = 4:7,
    if ~isempty(Tree{j,2}.ExtBasis),
        subplot(2,2,j-3);
        plot(Points,Tree{j,2}.ExtBasis(:,1));hold on;
        title(sprintf('Diffusion wavelet at scale %d',j));
    end;
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


function [vW,vDInvSqrt] = GraphNormalize(cW,cPrecision);

% Allocate memory
vW = cW;

% Compute the row sum
lD = sum( cW,2 );

lDInvSqrt = zeros(1,length(lD));
% Compute 'inverse' sqrt of lD
for lk = 1:length(lD)
    if abs(lD(lk)) > cPrecision,
        lDInvSqrt(lk) = 1/sqrt(lD(lk));
    else
        lDInvSqrt(lk) = 0;
    end;
end;

lDInvSqrt = diag(lDInvSqrt);

vW = lDInvSqrt*vW*lDInvSqrt;

%for lk = 1:size(vW,1)
%    vW(lk,lk) = 1-vW(lk,lk);
%end;

vDInvSqrt = lDInvSqrt;

return;