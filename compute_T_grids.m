% Compute the DW tree for the grid 
% Run make_T_grid.py before

clear all
close all
clc

T_grids = {'T_grid_16',	'T_grid_4',  'T_grid_8',
           'T_grid_sl_16',	'T_grid_sl_4',  'T_grid_sl_8'};
       
T_grids = {'T_grid_sl_16'};

GSOptions = struct('StopDensity',1,'Threshold',1e-2);
%opts = struct('Wavelets', true, 'WaveletPackets', true);%, 'OpThreshold', 1e-4, 'GSOptions', GSOptions);
opts = struct('Wavelets', true, 'OpThreshold', 1e-8, 'GSOptions', GSOptions);

for g = 1:length(T_grids),
    fn = T_grids{g};
    load(fn)
    Tree = DWPTree(T, 8, 1e-8, opts);

    mat_fn = ['Tree_', fn(3:end)]; 
    save(mat_fn, 'Tree')
end