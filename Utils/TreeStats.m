function Stats = TreeStats( Tree )

% Go through scales
for j = 1:length(Tree),    
    Stats.TSparsity(j)      = nnz(Tree{j,1}.T{1})/numel(Tree{j,1}.T{1});
    Stats.OpSparsity(j)     = nnz(Tree{j,1}.Op)/numel(Tree{j,1}.Op);
    Stats.BasisSparsity(j)  = nnz(Tree{j,1}.Basis)/numel(Tree{j,1}.Basis);
    Stats.Tnnz(j)           = nnz(Tree{j,1}.T{1});
    Stats.Opnnz(j)          = nnz(Tree{j,1}.Op);
    Stats.Basisnnz(j)       = nnz(Tree{j,1}.Basis);
end;

return;