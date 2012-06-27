function d = MatrixDensity( A, Thres )

if (nargin>1) && (Thres>0),
    A = abs(A)>Thres;
end;

d = nnz(A)/numel(A);

return;