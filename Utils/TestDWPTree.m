function TestDWPTree( T, Tree, Opts )

J = size(Tree,1);

% Test the accuracy of approximation of powers of the operator T
for j = 1:J,
    Stats.L2Op(j)     = max(svds(Tree{j,1}.ExtBasis*Tree{j,1}.T{1}*Tree{j,1}.ExtBasis'-T^(2^(j-1))));
    Stats.LinftyOp(j) = max(max(Tree{j,1}.ExtBasis*Tree{j,1}.T{1}*Tree{j,1}.ExtBasis'-T^(2^(j-1))));
    lS = svd(full(Tree{j,1}.Basis),0);
    Stats.BasisOrth(j)= max(lS)/min(lS);
end;

figure;semilogy(1:J,Stats.L2Op);title('L^2 approx of T^{2^j}');
figure;semilogy(1:J,Stats.LinftyOp);title('L^\infty approx of T^{2^j}');
figure;semilogy(1:J,Stats.BasisOrth-1);title('log(|Condition number of basis-1|)');

return;

