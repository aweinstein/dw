%% Construction
X=[0:127]/128;
A.vectors=X;
A.theta=.5;
A.maxdescend=int32(8);
B=covertree(A)

%% Find nearest example
C.vectors = X;
C.k       = int32(5);

[idxs,dists] = findnearest(X,B,C);

return;

%%
Y=X(:,B.levels(1:B.level_offsets(5))+1)

C.Y=Y;
C.within=.125;
C.level=int32(6);

Z=findwithin(X,B,C);

X(Z{2}+1)
