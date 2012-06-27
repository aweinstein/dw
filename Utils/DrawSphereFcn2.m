function h= DrawSphereFcn(Points, Fcn)
% function h = DrawSphereFcn(Points, Fcn)
%
% Draw a picture of the function Fcn on a sphere.
%
% In:
%    Points = Nx3 matrix of the N points comprising the sphere
%    Fcn    = Column vector given a function on those points.
% Out:
%    h      = handle to the generated surface
%
% Version History:
%    JCB 07/20/04
%

N = size(Points, 1);

% get the polargrid
[PX, PY, PZ] = sphere(200);
W = griddata3(Points(:,1), Points(:,2), Points(:,3), Fcn, PX, PY, PZ, 'nearest');

%figure;

colorbar
% C=gray(1024);
% C=C(100:900,:);
% colormap(C);
h=surf(PX, PY, PZ, W);
shading flat;
colorbar;

