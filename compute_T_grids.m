% Compute the DW tree for the grid 

T_grids = {'T_grid_16',	'T_grid_4',  'T_grid_8'};

GSOptions = struct('StopDensity',1,'Threshold',1e-2);
opts = struct('Wavelets', true, 'OpThreshold', 1e-4, 'GSOptions', GSOptions);

for g = 1:length(T_grids),
    fn = T_grids{g};
    load(fn)
    Tree = DWPTree(T, 15, 1e-4, GSOptions);
    mat_fn = ['Tree_', fn(3:end)]; 
    save(mat_fn, 'Tree')
end