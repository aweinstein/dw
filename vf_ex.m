clear all
close all
clc

load vf_ex1
fn = 'Tree_grid_ex1';

GSOptions = struct('StopDensity',1,'Threshold',1e-2);
opts = struct('Wavelets', true, 'OpThreshold', 1e-2, 'GSOptions', GSOptions);
tic
Tree = DWPTree(T, 15, 1e-10, GSOptions); 
toc
save(fn, 'Tree')
fprintf('Tree saved in %s\n', fn)