function [Y Z]=getlandmarks(X,B,l,w,r)

%X is D by N matrix of vectors; X(:,j) is j-th vector
%B what is returned from covertree
%l is as in C_l (int32); l is the level as in the cover tree notes
%w is how many beyond C_l we go (int32)
%r is distance
%Z is the cell array returned from findwithin; see below


N=B.params(1);
%number_duplicates=B.params(2);
%root=B.params(3);
%numlevels=B.params(4);
minlevel=B.params(5);
%maxlevel=B.params(6);
%maxnumlevels=B.params(7);

Y=X(:,B.levels(1:B.level_offsets(l-minlevel))+1);
C.Y=Y;
C.dist=r;
C.depth=int32(l+w);
C.space=[0:2*N-1];

Z=findwithin(X,B,C);

return;

