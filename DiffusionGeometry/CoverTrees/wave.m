X=[0:127]/128;
A.vectors=X;
A.theta=.5;
A.maxlevel=int32(8);
B=covertree(A)


Y=X(:,B.levels(1:B.level_offsets(5))+1)

C.Y=Y;
C.within=.125;
C.level=int32(6);

Z=findwithin(X,B,C);

X(Z{2}+1)
