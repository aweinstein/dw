epsilon = 0.001;

A = [epsilon,1-epsilon,0;1-epsilon,0,epsilon;0,0,1];

[lU,lS,lV] = svd(A);lS=diag(lS);

T = A/lS(1);

Tree = DWPTree (T, 20, 1e-8,struct('ScalePower',2,'StopFcns',1,'Wavelets',0));

T_highpower = T^10000;

[lUp,lSp,lVp]= svd(T_highpower);

figure;plot(Tree{size(Tree,1),1}.ExtBasis(:,1));hold on;plot(-lVp(:,1)','r');plot(T_highpower,'g');